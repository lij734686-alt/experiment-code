#define _CRT_SECURE_NO_WARNINGS

#include <pbc.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using Clock = std::chrono::high_resolution_clock;

static const char* DEFAULT_PAIRING_PARAM =
"type a\n"
"q 87807107996633125224377819847540498158068831994142082"
"1102865339926647563088022295707862517942266222142315585"
"8769582317459277713367317481324925129998224791\n"
"h 120160122648911460793888213667405342048029544012513118229196151310472072893597045311\n"
"r 730750818665451621361119245571504901405976559617\n"
"exp2 159\n"
"exp1 107\n"
"sign1 1\n"
"sign0 1\n";

struct Params {
    int N;      // authority number
    int ls;     // signature LSSS rows
    int le;     // decryption LSSS rows
    int attrs;  // attribute number
};

struct LSSS {
    element_t** M;
    int* rho;
    int rows;
    int cols;
};

struct System {
    pairing_ptr pairing;
    Params p;

    element_t g, K0, T0, beta1, beta2;
    element_t mu0;
    element_t* mu;

    element_t* alpha;
    element_t* Y;
    element_t** K;
    element_t** T;
};

struct SignKey {
    element_t* S0;
    element_t** S1;
    element_t** S2;
    element_t*** S3;
};

struct DecKey {
    element_t* D0;
    element_t** D1;
    element_t** D2;
    element_t*** D3;
};

struct Ciphertext {
    element_t C1;
    element_t C2;
    element_t C3;
    element_t sigma1;
    element_t sigma2;
    element_t sigma3;
};

static double elapsed_ms(Clock::time_point s, Clock::time_point e) {
    return std::chrono::duration<double, std::milli>(e - s).count();
}

static void die(void* p) {
    if (!p) {
        std::fprintf(stderr, "memory allocation failed\n");
        std::exit(1);
    }
}

static element_t* alloc_elements(int n) {
    element_t* x = static_cast<element_t*>(std::calloc((size_t)n, sizeof(element_t)));
    die(x);
    return x;
}

static element_t** alloc_matrix(int r, int c) {
    element_t** m = static_cast<element_t**>(std::calloc((size_t)r, sizeof(element_t*)));
    die(m);
    for (int i = 0; i < r; i++) m[i] = alloc_elements(c);
    return m;
}

static void init_lsss(LSSS* lsss, pairing_ptr pairing, int rows, int cols, int attrs) {
    lsss->rows = rows;
    lsss->cols = cols;
    lsss->rho = static_cast<int*>(std::calloc((size_t)rows, sizeof(int)));
    lsss->M = alloc_matrix(rows, cols);
    die(lsss->rho);

    for (int i = 0; i < rows; i++) {
        lsss->rho[i] = i % attrs;
        for (int j = 0; j < cols; j++) {
            element_init_Zr(lsss->M[i][j], pairing);
            element_set0(lsss->M[i][j]);
        }
    }

    // 简单可满足的 LSSS：第一行是 (1,0,...,0)
    element_set1(lsss->M[0][0]);
    for (int i = 1; i < rows; i++) {
        element_set1(lsss->M[i][i % cols]);
    }
}

static void dot_row(element_t out, pairing_ptr pairing, LSSS* lsss, int row, element_t* v) {
    element_set0(out);
    for (int j = 0; j < lsss->cols; j++) {
        element_t tmp;
        element_init_Zr(tmp, pairing);
        element_mul(tmp, lsss->M[row][j], v[j]);
        element_add(out, out, tmp);
        element_clear(tmp);
    }
}

static void init_system(System* sys, pairing_ptr pairing, Params p) {
    sys->pairing = pairing;
    sys->p = p;

    element_init_G1(sys->g, pairing);
    element_init_G1(sys->K0, pairing);
    element_init_G1(sys->T0, pairing);
    element_init_G1(sys->beta1, pairing);
    element_init_G1(sys->beta2, pairing);
    element_init_G1(sys->mu0, pairing);

    element_random(sys->g);
    element_random(sys->K0);
    element_random(sys->T0);
    element_random(sys->beta1);
    element_random(sys->beta2);
    element_random(sys->mu0);

    sys->mu = alloc_elements(p.attrs);
    for (int i = 0; i < p.attrs; i++) {
        element_init_G1(sys->mu[i], pairing);
        element_random(sys->mu[i]);
    }

    sys->alpha = alloc_elements(p.N);
    sys->Y = alloc_elements(p.N);
    sys->K = static_cast<element_t**>(std::calloc((size_t)p.N, sizeof(element_t*)));
    sys->T = static_cast<element_t**>(std::calloc((size_t)p.N, sizeof(element_t*)));
    die(sys->K);
    die(sys->T);

    element_t egg;
    element_init_GT(egg, pairing);
    pairing_apply(egg, sys->g, sys->g, pairing);

    for (int j = 0; j < p.N; j++) {
        element_init_Zr(sys->alpha[j], pairing);
        element_random(sys->alpha[j]);

        element_init_GT(sys->Y[j], pairing);
        element_pow_zn(sys->Y[j], egg, sys->alpha[j]);

        sys->K[j] = alloc_elements(p.attrs);
        sys->T[j] = alloc_elements(p.attrs);

        for (int i = 0; i < p.attrs; i++) {
            element_init_G1(sys->K[j][i], pairing);
            element_init_G1(sys->T[j][i], pairing);
            element_random(sys->K[j][i]);
            element_random(sys->T[j][i]);
        }
    }

    element_clear(egg);
}

static double sign_keygen(System* sys, LSSS* S, SignKey* sk) {
    auto start = Clock::now();

    Params p = sys->p;
    pairing_ptr pairing = sys->pairing;

    sk->S0 = alloc_elements(p.N);
    sk->S1 = static_cast<element_t**>(std::calloc((size_t)p.N, sizeof(element_t*)));
    sk->S2 = static_cast<element_t**>(std::calloc((size_t)p.N, sizeof(element_t*)));
    sk->S3 = static_cast<element_t***>(std::calloc((size_t)p.N, sizeof(element_t**)));
    die(sk->S1);
    die(sk->S2);
    die(sk->S3);

    for (int j = 0; j < p.N; j++) {
        element_t alpha_s1, alpha_s2;
        element_init_Zr(alpha_s1, pairing);
        element_init_Zr(alpha_s2, pairing);
        element_random(alpha_s1);
        element_sub(alpha_s2, sys->alpha[j], alpha_s1);

        element_init_G1(sk->S0[j], pairing);
        element_pow_zn(sk->S0[j], sys->g, alpha_s2);

        element_t* v = alloc_elements(S->cols);
        for (int k = 0; k < S->cols; k++) {
            element_init_Zr(v[k], pairing);
            if (k == 0) element_set(v[k], alpha_s1);
            else element_random(v[k]);
        }

        sk->S1[j] = alloc_elements(S->rows);
        sk->S2[j] = alloc_elements(S->rows);
        sk->S3[j] = static_cast<element_t**>(std::calloc((size_t)S->rows, sizeof(element_t*)));
        die(sk->S3[j]);

        for (int i = 0; i < S->rows; i++) {
            element_t lambda, r, part1, part2;
            element_init_Zr(lambda, pairing);
            element_init_Zr(r, pairing);
            element_init_G1(part1, pairing);
            element_init_G1(part2, pairing);
            element_random(r);

            dot_row(lambda, pairing, S, i, v);

            element_init_G1(sk->S1[j][i], pairing);
            element_init_G1(sk->S2[j][i], pairing);

            // S'_{j,i}=g^{lambda} K0 K_{j,rho(i)}^{r}
            element_pow_zn(part1, sys->g, lambda);
            element_pow_zn(part2, sys->K[j][S->rho[i]], r);
            element_mul(sk->S1[j][i], part1, sys->K0);
            element_mul(sk->S1[j][i], sk->S1[j][i], part2);

            // S''_{j,i}=g^r
            element_pow_zn(sk->S2[j][i], sys->g, r);

            // S'''_{j,i,x}=K_{j,x}^r
            sk->S3[j][i] = alloc_elements(p.attrs);
            for (int x = 0; x < p.attrs; x++) {
                element_init_G1(sk->S3[j][i][x], pairing);
                if (x == S->rho[i]) element_set1(sk->S3[j][i][x]);
                else element_pow_zn(sk->S3[j][i][x], sys->K[j][x], r);
            }

            element_clear(lambda);
            element_clear(r);
            element_clear(part1);
            element_clear(part2);
        }

        for (int k = 0; k < S->cols; k++) element_clear(v[k]);
        std::free(v);
        element_clear(alpha_s1);
        element_clear(alpha_s2);
    }

    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static double dec_keygen(System* sys, LSSS* D, DecKey* dk) {
    auto start = Clock::now();

    Params p = sys->p;
    pairing_ptr pairing = sys->pairing;

    dk->D0 = alloc_elements(p.N);
    dk->D1 = static_cast<element_t**>(std::calloc((size_t)p.N, sizeof(element_t*)));
    dk->D2 = static_cast<element_t**>(std::calloc((size_t)p.N, sizeof(element_t*)));
    dk->D3 = static_cast<element_t***>(std::calloc((size_t)p.N, sizeof(element_t**)));
    die(dk->D1);
    die(dk->D2);
    die(dk->D3);

    for (int j = 0; j < p.N; j++) {
        element_t alpha_e1, alpha_e2;
        element_init_Zr(alpha_e1, pairing);
        element_init_Zr(alpha_e2, pairing);
        element_random(alpha_e1);
        element_sub(alpha_e2, sys->alpha[j], alpha_e1);

        element_init_G1(dk->D0[j], pairing);
        element_pow_zn(dk->D0[j], sys->g, alpha_e2);

        element_t* v = alloc_elements(D->cols);
        for (int k = 0; k < D->cols; k++) {
            element_init_Zr(v[k], pairing);
            if (k == 0) element_set(v[k], alpha_e1);
            else element_random(v[k]);
        }

        dk->D1[j] = alloc_elements(D->rows);
        dk->D2[j] = alloc_elements(D->rows);
        dk->D3[j] = static_cast<element_t**>(std::calloc((size_t)D->rows, sizeof(element_t*)));
        die(dk->D3[j]);

        for (int i = 0; i < D->rows; i++) {
            element_t lambda, t, part1, part2;
            element_init_Zr(lambda, pairing);
            element_init_Zr(t, pairing);
            element_init_G1(part1, pairing);
            element_init_G1(part2, pairing);
            element_random(t);

            dot_row(lambda, pairing, D, i, v);

            element_init_G1(dk->D1[j][i], pairing);
            element_init_G1(dk->D2[j][i], pairing);

            // D'_{j,i}=g^{lambda} T0 T_{j,rho(i)}^t
            element_pow_zn(part1, sys->g, lambda);
            element_pow_zn(part2, sys->T[j][D->rho[i]], t);
            element_mul(dk->D1[j][i], part1, sys->T0);
            element_mul(dk->D1[j][i], dk->D1[j][i], part2);

            // D''_{j,i}=g^t
            element_pow_zn(dk->D2[j][i], sys->g, t);

            // D'''_{j,i,x}=T_{j,x}^t
            dk->D3[j][i] = alloc_elements(p.attrs);
            for (int x = 0; x < p.attrs; x++) {
                element_init_G1(dk->D3[j][i][x], pairing);
                if (x == D->rho[i]) element_set1(dk->D3[j][i][x]);
                else element_pow_zn(dk->D3[j][i][x], sys->T[j][x], t);
            }

            element_clear(lambda);
            element_clear(t);
            element_clear(part1);
            element_clear(part2);
        }

        for (int k = 0; k < D->cols; k++) element_clear(v[k]);
        std::free(v);
        element_clear(alpha_e1);
        element_clear(alpha_e2);
    }

    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static double encryption(System* sys, Ciphertext* ct) {
    auto start = Clock::now();

    Params p = sys->p;
    pairing_ptr pairing = sys->pairing;

    element_t beta, tempG;
    element_init_Zr(beta, pairing);
    element_init_G1(tempG, pairing);
    element_random(beta);

    element_init_G1(ct->C1, pairing);
    element_init_G1(ct->C2, pairing);
    element_init_GT(ct->C3, pairing);

    // C1=g^beta
    element_pow_zn(ct->C1, sys->g, beta);

    // C2=(T0 prod T_{j,i})^beta
    element_set(tempG, sys->T0);
    for (int j = 0; j < p.N; j++) {
        element_mul(tempG, tempG, sys->T[j][0]);
    }
    element_pow_zn(ct->C2, tempG, beta);

    // C3=H2((prod Y_j)^beta, sigma1, eta) xor M，这里构造 GT 部分
    element_t prodY;
    element_init_GT(prodY, pairing);
    element_set1(prodY);
    for (int j = 0; j < p.N; j++) {
        element_mul(prodY, prodY, sys->Y[j]);
    }
    element_pow_zn(ct->C3, prodY, beta);

    element_clear(beta);
    element_clear(tempG);
    element_clear(prodY);

    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static double signature(System* sys, LSSS* S, SignKey* sk, Ciphertext* ct) {
    auto start = Clock::now();

    Params p = sys->p;
    pairing_ptr pairing = sys->pairing;

    element_t beta, gamma, nu, temp, partial1, partial2;
    element_init_Zr(beta, pairing);
    element_init_Zr(gamma, pairing);
    element_init_Zr(nu, pairing);
    element_init_G1(temp, pairing);
    element_init_G1(partial1, pairing);
    element_init_G1(partial2, pairing);

    element_random(beta);
    element_random(gamma);
    element_random(nu);

    element_init_G1(ct->sigma1, pairing);
    element_init_G1(ct->sigma2, pairing);
    element_init_G1(ct->sigma3, pairing);

    // ES partial signature sigma1'=g^nu prod S''_{j,0}
    element_pow_zn(partial1, sys->g, nu);
    for (int j = 0; j < p.N; j++) {
        element_mul(partial1, partial1, sk->S2[j][0]);
    }

    // ES partial signature sigma2'=K0 prod K_{j,rho(0)} prod S'_{j,0}
    element_set(partial2, sys->K0);
    for (int j = 0; j < p.N; j++) {
        element_mul(partial2, partial2, sys->K[j][S->rho[0]]);
        element_mul(partial2, partial2, sk->S1[j][0]);
    }

    // sigma1=g^{beta*gamma}
    element_mul(gamma, gamma, beta);
    element_pow_zn(ct->sigma1, sys->g, gamma);

    // sigma2=sigma1'
    element_set(ct->sigma2, partial1);

    // sigma3=beta1^delta beta2^gamma mu0 prod mu_k^{m_k} sigma2' prod S0_j
    element_random(ct->sigma3);
    element_mul(ct->sigma3, ct->sigma3, partial2);
    for (int j = 0; j < p.N; j++) {
        element_mul(ct->sigma3, ct->sigma3, sk->S0[j]);
    }

    element_clear(beta);
    element_clear(gamma);
    element_clear(nu);
    element_clear(temp);
    element_clear(partial1);
    element_clear(partial2);

    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static double verification(System* sys, Ciphertext* ct) {
    auto start = Clock::now();

    pairing_ptr pairing = sys->pairing;

    element_t V1, V2, tempG;
    element_init_GT(V1, pairing);
    element_init_GT(V2, pairing);
    element_init_G1(tempG, pairing);

    // Verify Eq.(7): local verification mainly computes pairings.
    element_random(tempG);
    pairing_apply(V1, tempG, ct->sigma1, pairing);

    element_random(tempG);
    pairing_apply(V2, tempG, ct->C1, pairing);

    pairing_apply(V1, sys->K0, ct->sigma2, pairing);
    pairing_apply(V2, ct->sigma3, sys->g, pairing);

    element_clear(V1);
    element_clear(V2);
    element_clear(tempG);

    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static double decryption(System* sys, LSSS* D, DecKey* dk, Ciphertext* ct) {
    auto start = Clock::now();

    Params p = sys->p;
    pairing_ptr pairing = sys->pairing;

    element_t C1p, C2p, tempG, tempGT, Z;
    element_init_G1(C1p, pairing);
    element_init_GT(C2p, pairing);
    element_init_G1(tempG, pairing);
    element_init_GT(tempGT, pairing);
    element_init_GT(Z, pairing);

    element_set1(C1p);
    element_set1(C2p);

    // ES computes C1' = prod D'_{j,0} prod D'''_{j,0,x}
    for (int j = 0; j < p.N; j++) {
        element_mul(C1p, C1p, dk->D1[j][0]);
        for (int x = 0; x < p.attrs; x++) {
            if (x != D->rho[0]) {
                element_mul(C1p, C1p, dk->D3[j][0][x]);
            }
        }
    }

    // ES computes C2' = prod e(C2,D''_{j,0})
    for (int j = 0; j < p.N; j++) {
        pairing_apply(tempGT, ct->C2, dk->D2[j][0], pairing);
        element_mul(C2p, C2p, tempGT);
    }

    // DU computes e(C1, C1' prod D0_j) / C2'
    for (int j = 0; j < p.N; j++) {
        element_mul(C1p, C1p, dk->D0[j]);
    }
    pairing_apply(Z, ct->C1, C1p, pairing);
    element_div(Z, Z, C2p);

    element_clear(C1p);
    element_clear(C2p);
    element_clear(tempG);
    element_clear(tempGT);
    element_clear(Z);

    auto end = Clock::now();
    return elapsed_ms(start, end);
}

int main(int argc, char** argv) {
    Params p = { 5, 10, 10, 10 };

    if (argc > 1) p.N = std::atoi(argv[1]);
    if (argc > 2) p.ls = std::atoi(argv[2]);
    if (argc > 3) p.le = std::atoi(argv[3]);
    if (argc > 4) p.attrs = std::atoi(argv[4]);

    pairing_t pairing;
    pairing_init_set_buf(pairing, DEFAULT_PAIRING_PARAM, std::strlen(DEFAULT_PAIRING_PARAM));

    std::printf("blocks,dec_keygen_ms,sign_keygen_ms,encryption_ms,signature_ms,decryption_ms,verification_ms\n");

    for (int blocks = 50; blocks <= 1000; blocks += 50) {
        System sys;
        LSSS S, D;
        SignKey sk;
        DecKey dk;

        init_system(&sys, pairing, p);
        init_lsss(&S, pairing, p.ls, p.attrs, p.attrs);
        init_lsss(&D, pairing, p.le, p.attrs, p.attrs);

        // 用户数量固定为 1，所以密钥只生成一次
        double t_sign_keygen = sign_keygen(&sys, &S, &sk);
        double t_dec_keygen = dec_keygen(&sys, &D, &dk);

        double t_encryption = 0.0;
        double t_signature = 0.0;
        double t_decryption = 0.0;
        double t_verification = 0.0;

        for (int b = 0; b < blocks; b++) {
            Ciphertext ct;

            t_encryption += encryption(&sys, &ct);
            t_signature += signature(&sys, &S, &sk, &ct);
            t_verification += verification(&sys, &ct);
            t_decryption += decryption(&sys, &D, &dk, &ct);
        }

        std::printf("%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            blocks,
            t_dec_keygen,
            t_sign_keygen,
            t_encryption,
            t_signature,
            t_decryption,
            t_verification);
    }

    pairing_clear(pairing);
    return 0;
}