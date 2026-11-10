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

struct SchemeParams {
    int n_attrs;       // 属性数量 n
    int l_rows;        // LSSS 矩阵行数 l
    int poly_degree;   // 多项式阶数 L
    int file_count;    // 文件数量 N
};

struct LsssPolicy {
    element_t** M;
    int* rho;
    int rows;
    int cols;
};

struct PublicContext {
    pairing_ptr pairing;
    SchemeParams params;

    element_t g, u;
    element_t alpha_s, alpha_A, beta;
    element_t g_beta;
    element_t PK_T;
    element_t h0;
    element_t* h;
    element_t* P;
};

struct AuditorKey {
    element_t K_A;
    element_t L_A;
    element_t* K_A_x;
};

struct FileCiphertext {
    element_t x_j, X_j;
    element_t theta_j, Y_j;
    element_t HF_j, h0_j;
    element_t s_j;
    element_t C_hat_prime;
    element_t* C_i;
    element_t* D_i;
};

struct SigningKey {
    element_t* K_i;
    element_t* K_i_prime;
    element_t** K_i_x_double_prime;
};

struct BlockSignature {
    element_t sigma0;
    element_t sigma1;
    element_t sigma2;
};

static double elapsed_ms(Clock::time_point s, Clock::time_point e) {
    return std::chrono::duration<double, std::milli>(e - s).count();
}

static void die_alloc(void* p) {
    if (!p) {
        std::fprintf(stderr, "memory allocation failed\n");
        std::exit(1);
    }
}

static element_t* alloc_elements(int count) {
    element_t* p = static_cast<element_t*>(std::calloc((size_t)count, sizeof(element_t)));
    die_alloc(p);
    return p;
}

static void init_public_context(PublicContext* ctx, pairing_ptr pairing, SchemeParams params) {
    ctx->pairing = pairing;
    ctx->params = params;

    element_init_G1(ctx->g, pairing);
    element_init_G1(ctx->u, pairing);
    element_init_Zr(ctx->alpha_s, pairing);
    element_init_Zr(ctx->alpha_A, pairing);
    element_init_Zr(ctx->beta, pairing);
    element_init_G1(ctx->g_beta, pairing);
    element_init_GT(ctx->PK_T, pairing);
    element_init_G1(ctx->h0, pairing);

    element_random(ctx->g);
    element_random(ctx->u);
    element_random(ctx->alpha_s);
    element_random(ctx->alpha_A);
    element_random(ctx->beta);
    element_random(ctx->h0);

    element_pow_zn(ctx->g_beta, ctx->g, ctx->beta);

    element_t egg;
    element_init_GT(egg, pairing);
    pairing_apply(egg, ctx->g, ctx->g, pairing);
    element_pow_zn(ctx->PK_T, egg, ctx->alpha_s);
    element_clear(egg);

    ctx->h = alloc_elements(params.n_attrs);
    for (int i = 0; i < params.n_attrs; i++) {
        element_init_G1(ctx->h[i], pairing);
        element_random(ctx->h[i]);
    }

    ctx->P = alloc_elements(params.poly_degree + 1);
    for (int i = 0; i <= params.poly_degree; i++) {
        element_t x;
        element_init_Zr(x, pairing);
        element_random(x);
        element_init_G1(ctx->P[i], pairing);
        element_pow_zn(ctx->P[i], ctx->g, x);
        element_clear(x);
    }
}

static void clear_public_context(PublicContext* ctx) {
    element_clear(ctx->g);
    element_clear(ctx->u);
    element_clear(ctx->alpha_s);
    element_clear(ctx->alpha_A);
    element_clear(ctx->beta);
    element_clear(ctx->g_beta);
    element_clear(ctx->PK_T);
    element_clear(ctx->h0);

    for (int i = 0; i < ctx->params.n_attrs; i++) element_clear(ctx->h[i]);
    std::free(ctx->h);

    for (int i = 0; i <= ctx->params.poly_degree; i++) element_clear(ctx->P[i]);
    std::free(ctx->P);
}

static void init_policy(LsssPolicy* policy, pairing_ptr pairing, int rows, int cols, int n_attrs) {
    policy->rows = rows;
    policy->cols = cols;

    policy->rho = static_cast<int*>(std::calloc((size_t)rows, sizeof(int)));
    policy->M = static_cast<element_t**>(std::calloc((size_t)rows, sizeof(element_t*)));
    die_alloc(policy->rho);
    die_alloc(policy->M);

    for (int i = 0; i < rows; i++) {
        policy->rho[i] = i % n_attrs;
        policy->M[i] = alloc_elements(cols);
        for (int j = 0; j < cols; j++) {
            element_init_Zr(policy->M[i][j], pairing);
            element_set0(policy->M[i][j]);
        }
    }

    // 构造简单可满足 LSSS 矩阵：第一行是 (1,0,...,0)
    element_set1(policy->M[0][0]);
    for (int i = 1; i < rows; i++) {
        element_set1(policy->M[i][i % cols]);
    }
}

static void clear_policy(LsssPolicy* policy) {
    for (int i = 0; i < policy->rows; i++) {
        for (int j = 0; j < policy->cols; j++) element_clear(policy->M[i][j]);
        std::free(policy->M[i]);
    }
    std::free(policy->M);
    std::free(policy->rho);
}

static void dot_row_vector(element_t out, pairing_ptr pairing, LsssPolicy* policy, int row, element_t* vec) {
    element_set0(out);
    for (int j = 0; j < policy->cols; j++) {
        element_t tmp;
        element_init_Zr(tmp, pairing);
        element_mul(tmp, policy->M[row][j], vec[j]);
        element_add(out, out, tmp);
        element_clear(tmp);
    }
}

static double dec_keygen(PublicContext* ctx, AuditorKey* ak) {
    auto start = Clock::now();

    int n = ctx->params.n_attrs;
    pairing_ptr pairing = ctx->pairing;

    element_t t_A, part1, part2;
    element_init_Zr(t_A, pairing);
    element_init_G1(part1, pairing);
    element_init_G1(part2, pairing);
    element_random(t_A);

    element_init_G1(ak->K_A, pairing);
    element_init_G1(ak->L_A, pairing);
    ak->K_A_x = alloc_elements(n);

    // K_A = g^{alpha_A}(g^beta)^{t_A}
    element_pow_zn(part1, ctx->g, ctx->alpha_A);
    element_pow_zn(part2, ctx->g_beta, t_A);
    element_mul(ak->K_A, part1, part2);

    // L_A = g^{t_A}
    element_pow_zn(ak->L_A, ctx->g, t_A);

    // K_{A,x}=h_x^{t_A}
    for (int x = 0; x < n; x++) {
        element_init_G1(ak->K_A_x[x], pairing);
        element_pow_zn(ak->K_A_x[x], ctx->h[x], t_A);
    }

    element_clear(t_A);
    element_clear(part1);
    element_clear(part2);

    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static void clear_auditor_key(AuditorKey* ak, int n) {
    element_clear(ak->K_A);
    element_clear(ak->L_A);
    for (int i = 0; i < n; i++) element_clear(ak->K_A_x[i]);
    std::free(ak->K_A_x);
}

static double encryption(PublicContext* ctx, LsssPolicy* policy, FileCiphertext* ct) {
    auto start = Clock::now();

    int l = policy->rows;
    int cols = policy->cols;
    pairing_ptr pairing = ctx->pairing;

    element_init_Zr(ct->x_j, pairing);
    element_init_G1(ct->X_j, pairing);
    element_init_Zr(ct->theta_j, pairing);
    element_init_G1(ct->Y_j, pairing);
    element_init_G1(ct->HF_j, pairing);
    element_init_G1(ct->h0_j, pairing);
    element_init_Zr(ct->s_j, pairing);
    element_init_G1(ct->C_hat_prime, pairing);

    element_random(ct->x_j);
    element_pow_zn(ct->X_j, ctx->g, ct->x_j);

    element_random(ct->theta_j);
    element_pow_zn(ct->Y_j, ctx->g, ct->theta_j);

    element_random(ct->HF_j);
    element_random(ct->h0_j);
    element_random(ct->s_j);
    element_pow_zn(ct->C_hat_prime, ctx->g, ct->s_j);

    element_t* y = alloc_elements(cols);
    for (int i = 0; i < cols; i++) {
        element_init_Zr(y[i], pairing);
        if (i == 0) element_set(y[i], ct->s_j);
        else element_random(y[i]);
    }

    ct->C_i = alloc_elements(l);
    ct->D_i = alloc_elements(l);

    for (int i = 0; i < l; i++) {
        element_t lambda, r, left, right, inv;
        element_init_Zr(lambda, pairing);
        element_init_Zr(r, pairing);
        element_init_G1(left, pairing);
        element_init_G1(right, pairing);
        element_init_G1(inv, pairing);

        element_random(r);
        dot_row_vector(lambda, pairing, policy, i, y);

        // C_i=(g^beta)^lambda * h_{rho(i)}^{-r}
        element_init_G1(ct->C_i[i], pairing);
        element_init_G1(ct->D_i[i], pairing);

        element_pow_zn(left, ctx->g_beta, lambda);
        element_pow_zn(right, ctx->h[policy->rho[i]], r);
        element_invert(inv, right);
        element_mul(ct->C_i[i], left, inv);

        // D_i=g^r
        element_pow_zn(ct->D_i[i], ctx->g, r);

        element_clear(lambda);
        element_clear(r);
        element_clear(left);
        element_clear(right);
        element_clear(inv);
    }

    for (int i = 0; i < cols; i++) element_clear(y[i]);
    std::free(y);

    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static void clear_file_ciphertext(FileCiphertext* ct, LsssPolicy* policy) {
    element_clear(ct->x_j);
    element_clear(ct->X_j);
    element_clear(ct->theta_j);
    element_clear(ct->Y_j);
    element_clear(ct->HF_j);
    element_clear(ct->h0_j);
    element_clear(ct->s_j);
    element_clear(ct->C_hat_prime);

    for (int i = 0; i < policy->rows; i++) {
        element_clear(ct->C_i[i]);
        element_clear(ct->D_i[i]);
    }
    std::free(ct->C_i);
    std::free(ct->D_i);
}

static double sign_keygen(PublicContext* ctx, LsssPolicy* policy, FileCiphertext* ct, SigningKey* sk) {
    auto start = Clock::now();

    int l = policy->rows;
    int n = ctx->params.n_attrs;
    int cols = policy->cols;
    pairing_ptr pairing = ctx->pairing;

    element_t* v = alloc_elements(cols);
    for (int i = 0; i < cols; i++) {
        element_init_Zr(v[i], pairing);
        if (i == 0) element_set(v[i], ctx->alpha_s);
        else element_random(v[i]);
    }

    sk->K_i = alloc_elements(l);
    sk->K_i_prime = alloc_elements(l);
    sk->K_i_x_double_prime = static_cast<element_t**>(std::calloc((size_t)l, sizeof(element_t*)));
    die_alloc(sk->K_i_x_double_prime);

    for (int i = 0; i < l; i++) {
        element_t lambda, u_i, base, part1, part2;
        element_init_Zr(lambda, pairing);
        element_init_Zr(u_i, pairing);
        element_init_G1(base, pairing);
        element_init_G1(part1, pairing);
        element_init_G1(part2, pairing);

        element_random(u_i);
        dot_row_vector(lambda, pairing, policy, i, v);

        // K_i=g^lambda * (h0_j h_{rho(i)})^{u_i}
        element_init_G1(sk->K_i[i], pairing);
        element_mul(base, ct->h0_j, ctx->h[policy->rho[i]]);
        element_pow_zn(part1, ctx->g, lambda);
        element_pow_zn(part2, base, u_i);
        element_mul(sk->K_i[i], part1, part2);

        // K_i'=g^{u_i}
        element_init_G1(sk->K_i_prime[i], pairing);
        element_pow_zn(sk->K_i_prime[i], ctx->g, u_i);

        // K''_{i,x}=h_x^{u_i}
        sk->K_i_x_double_prime[i] = alloc_elements(n);
        for (int x = 0; x < n; x++) {
            element_init_G1(sk->K_i_x_double_prime[i][x], pairing);
            if (x == policy->rho[i]) element_set1(sk->K_i_x_double_prime[i][x]);
            else element_pow_zn(sk->K_i_x_double_prime[i][x], ctx->h[x], u_i);
        }

        element_clear(lambda);
        element_clear(u_i);
        element_clear(base);
        element_clear(part1);
        element_clear(part2);
    }

    for (int i = 0; i < cols; i++) element_clear(v[i]);
    std::free(v);

    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static void clear_signing_key(SigningKey* sk, LsssPolicy* policy, int n) {
    for (int i = 0; i < policy->rows; i++) {
        element_clear(sk->K_i[i]);
        element_clear(sk->K_i_prime[i]);
        for (int x = 0; x < n; x++) element_clear(sk->K_i_x_double_prime[i][x]);
        std::free(sk->K_i_x_double_prime[i]);
    }
    std::free(sk->K_i);
    std::free(sk->K_i_prime);
    std::free(sk->K_i_x_double_prime);
}

static void compute_R(element_t R, PublicContext* ctx) {
    element_set(R, ctx->P[0]);
    for (int i = 1; i <= ctx->params.poly_degree; i++) {
        element_t h, part;
        element_init_Zr(h, ctx->pairing);
        element_init_G1(part, ctx->pairing);
        element_random(h);
        element_pow_zn(part, ctx->P[i], h);
        element_mul(R, R, part);
        element_clear(h);
        element_clear(part);
    }
}

static double block_signature(PublicContext* ctx, LsssPolicy* policy, FileCiphertext* ct,
    SigningKey* sk, BlockSignature* sigs, int blocks) {
    auto start = Clock::now();

    int n = ctx->params.n_attrs;
    pairing_ptr pairing = ctx->pairing;

    for (int t = 0; t < blocks; t++) {
        element_t xi, r, R, HW, tmp, product, msg, rpart;

        element_init_Zr(xi, pairing);
        element_init_Zr(r, pairing);
        element_init_G1(R, pairing);
        element_init_G1(HW, pairing);
        element_init_G1(tmp, pairing);
        element_init_G1(product, pairing);
        element_init_G1(msg, pairing);
        element_init_G1(rpart, pairing);

        element_random(xi);
        element_random(r);
        compute_R(R, ctx);

        element_init_G1(sigs[t].sigma0, pairing);
        element_init_G1(sigs[t].sigma1, pairing);
        element_init_G1(sigs[t].sigma2, pairing);

        // sigma0=g^r
        element_pow_zn(sigs[t].sigma0, ctx->g, r);

        // 简化策略下 I={0}, omega=1
        element_pow_zn(tmp, ctx->g, xi);
        element_mul(sigs[t].sigma1, sk->K_i_prime[0], tmp);

        element_set(HW, ct->h0_j);
        for (int x = 0; x < n; x++) element_mul(HW, HW, ctx->h[x]);

        element_set(product, sk->K_i[0]);
        for (int x = 0; x < n; x++) {
            if (x != policy->rho[0]) {
                element_mul(product, product, sk->K_i_x_double_prime[0][x]);
            }
        }

        element_pow_zn(tmp, HW, xi);
        element_mul(product, product, tmp);
        element_mul(product, product, ct->HF_j);

        element_random(msg);
        element_pow_zn(msg, msg, ct->x_j);
        element_mul(product, product, msg);

        element_pow_zn(rpart, R, r);
        element_mul(sigs[t].sigma2, product, rpart);

        element_clear(xi);
        element_clear(r);
        element_clear(R);
        element_clear(HW);
        element_clear(tmp);
        element_clear(product);
        element_clear(msg);
        element_clear(rpart);
    }

    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static void clear_block_signatures(BlockSignature* sigs, int blocks) {
    for (int i = 0; i < blocks; i++) {
        element_clear(sigs[i].sigma0);
        element_clear(sigs[i].sigma1);
        element_clear(sigs[i].sigma2);
    }
}

static double decrypt_authorization(PublicContext* ctx, LsssPolicy* policy, AuditorKey* ak, FileCiphertext* ct) {
    auto start = Clock::now();

    pairing_ptr pairing = ctx->pairing;
    element_t numerator, denominator, p1, p2, tmp, Z;

    element_init_GT(numerator, pairing);
    element_init_GT(denominator, pairing);
    element_init_GT(p1, pairing);
    element_init_GT(p2, pairing);
    element_init_GT(tmp, pairing);
    element_init_GT(Z, pairing);

    element_set1(denominator);

    pairing_apply(numerator, ct->C_hat_prime, ak->K_A, pairing);

    for (int i = 0; i < policy->rows; i++) {
        pairing_apply(p1, ct->C_i[i], ak->L_A, pairing);
        pairing_apply(p2, ct->D_i[i], ak->K_A_x[policy->rho[i]], pairing);
        element_mul(tmp, p1, p2);
        element_mul(denominator, denominator, tmp);
    }

    element_invert(denominator, denominator);
    element_mul(Z, numerator, denominator);

    element_clear(numerator);
    element_clear(denominator);
    element_clear(p1);
    element_clear(p2);
    element_clear(tmp);
    element_clear(Z);

    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static double verify_audit(PublicContext* ctx, FileCiphertext* ct, BlockSignature* sigs, int blocks) {
    auto start = Clock::now();

    pairing_ptr pairing = ctx->pairing;
    element_t left, right, tmpGT, tmpG1, exp;

    element_init_GT(left, pairing);
    element_init_GT(right, pairing);
    element_init_GT(tmpGT, pairing);
    element_init_G1(tmpG1, pairing);
    element_init_Zr(exp, pairing);

    // 左边 e(delta2,g)，这里用随机聚合元素模拟 delta2
    element_random(tmpG1);
    pairing_apply(left, tmpG1, ctx->g, pairing);

    element_set1(right);

    element_random(exp);
    element_pow_zn(tmpGT, ctx->PK_T, exp);
    element_mul(right, right, tmpGT);

    // e(HW,delta1)
    element_t HW, delta1;
    element_init_G1(HW, pairing);
    element_init_G1(delta1, pairing);
    element_set(HW, ct->h0_j);
    for (int x = 0; x < ctx->params.n_attrs; x++) element_mul(HW, HW, ctx->h[x]);
    element_random(delta1);
    pairing_apply(tmpGT, HW, delta1, pairing);
    element_mul(right, right, tmpGT);

    // e(HF^V,g)
    element_random(exp);
    element_pow_zn(tmpG1, ct->HF_j, exp);
    pairing_apply(tmpGT, tmpG1, ctx->g, pairing);
    element_mul(right, right, tmpGT);

    // e(product, X_j)
    element_random(tmpG1);
    pairing_apply(tmpGT, tmpG1, ct->X_j, pairing);
    element_mul(right, right, tmpGT);

    // Omega_R = prod e(R, sigma0)
    for (int t = 0; t < blocks; t++) {
        element_t R;
        element_init_G1(R, pairing);
        compute_R(R, ctx);
        pairing_apply(tmpGT, R, sigs[t].sigma0, pairing);
        element_mul(right, right, tmpGT);
        element_clear(R);
    }

    element_clear(left);
    element_clear(right);
    element_clear(tmpGT);
    element_clear(tmpG1);
    element_clear(exp);
    element_clear(HW);
    element_clear(delta1);

    auto end = Clock::now();
    return elapsed_ms(start, end);
}

int main(int argc, char** argv) {
    SchemeParams params = { 10, 10, 5, 1 };

    if (argc > 1) params.n_attrs = std::atoi(argv[1]);
    if (argc > 2) params.l_rows = std::atoi(argv[2]);
    if (argc > 3) params.poly_degree = std::atoi(argv[3]);
    if (argc > 4) params.file_count = std::atoi(argv[4]);

    pairing_t pairing;
    pairing_init_set_buf(pairing, DEFAULT_PAIRING_PARAM, std::strlen(DEFAULT_PAIRING_PARAM));

    int block_sizes[] = { 50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000 };
    int block_size_count = sizeof(block_sizes) / sizeof(block_sizes[0]);

    std::printf("blocks,dec_keygen_ms,sign_keygen_ms,encryption_ms,signature_ms,decryption_ms,verification_ms\n");

    for (int k = 0; k < block_size_count; k++) {
        int blocks = block_sizes[k];

        PublicContext ctx;
        LsssPolicy policy;
        AuditorKey ak;
        FileCiphertext ct;
        SigningKey sk;

        init_public_context(&ctx, pairing, params);
        init_policy(&policy, pairing, params.l_rows, params.n_attrs, params.n_attrs);

        BlockSignature* sigs = static_cast<BlockSignature*>(std::calloc((size_t)blocks, sizeof(BlockSignature)));
        die_alloc(sigs);

        double t_dec_keygen = dec_keygen(&ctx, &ak);
        double t_encryption = encryption(&ctx, &policy, &ct);
        double t_sign_keygen = sign_keygen(&ctx, &policy, &ct, &sk);
        double t_signature = block_signature(&ctx, &policy, &ct, &sk, sigs, blocks);
        double t_decryption = decrypt_authorization(&ctx, &policy, &ak, &ct);
        double t_verification = verify_audit(&ctx, &ct, sigs, blocks);

        std::printf("%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            blocks,
            t_dec_keygen,
            t_sign_keygen,
            t_encryption,
            t_signature,
            t_decryption,
            t_verification);

        clear_block_signatures(sigs, blocks);
        std::free(sigs);
        clear_signing_key(&sk, &policy, params.n_attrs);
        clear_file_ciphertext(&ct, &policy);
        clear_auditor_key(&ak, params.n_attrs);
        clear_policy(&policy);
        clear_public_context(&ctx);
    }

    pairing_clear(pairing);
    return 0;
}