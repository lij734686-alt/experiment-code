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
    int n;   // 文档关键词数量
    int m;   // 查询关键词数量
};

struct System {
    pairing_ptr pairing;

    element_t g, g1, g2;
    element_t Y, EK;

    element_t b1, b2, a1, x1;
    element_t alpha, lambda;

    element_t skcs, pkcs;
    element_t skts, pkts;
    element_t skp, pkp;
};

struct DoctorKey {
    element_t D1, D2, D3, D4, D5;
    element_t pkd;
};

struct Ciphertext {
    element_t* C1;   // C1[0], C1[1..n]
    element_t C2;
    element_t C3;
};

struct Trapdoor {
    element_t* T1;   // T1[0], T1[1..n]
    element_t T2;
    element_t T3;
    element_t T4;
    element_t T5;
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

static element_t* alloc_elements(int count) {
    element_t* arr = static_cast<element_t*>(std::calloc((size_t)count, sizeof(element_t)));
    die(arr);
    return arr;
}

static void poly_f(element_t out, System* sys, element_t x) {
    element_t tmp;
    element_init_Zr(tmp, sys->pairing);
    element_mul(tmp, sys->a1, x);
    element_add(out, sys->b1, tmp);
    element_clear(tmp);
}

static void init_system(System* sys, pairing_ptr pairing) {
    sys->pairing = pairing;

    element_init_G1(sys->g, pairing);
    element_init_G1(sys->g1, pairing);
    element_init_G1(sys->g2, pairing);
    element_init_GT(sys->Y, pairing);
    element_init_G1(sys->EK, pairing);

    element_init_Zr(sys->b1, pairing);
    element_init_Zr(sys->b2, pairing);
    element_init_Zr(sys->a1, pairing);
    element_init_Zr(sys->x1, pairing);
    element_init_Zr(sys->alpha, pairing);
    element_init_Zr(sys->lambda, pairing);

    element_random(sys->g);
    element_random(sys->b1);
    element_random(sys->b2);
    element_random(sys->a1);
    element_random(sys->x1);
    element_random(sys->alpha);
    element_random(sys->lambda);

    // g1=g^{b1}, g2=g^{b2}
    element_pow_zn(sys->g1, sys->g, sys->b1);
    element_pow_zn(sys->g2, sys->g, sys->b2);

    // Y=e(g,g)^alpha
    element_t egg;
    element_init_GT(egg, pairing);
    pairing_apply(egg, sys->g, sys->g, pairing);
    element_pow_zn(sys->Y, egg, sys->alpha);

    // EK=g^{f(x1)/b1}
    element_t fx1, inv_b1, exp;
    element_init_Zr(fx1, pairing);
    element_init_Zr(inv_b1, pairing);
    element_init_Zr(exp, pairing);

    poly_f(fx1, sys, sys->x1);
    element_invert(inv_b1, sys->b1);
    element_mul(exp, fx1, inv_b1);
    element_pow_zn(sys->EK, sys->g, exp);

    element_clear(egg);
    element_clear(fx1);
    element_clear(inv_b1);
    element_clear(exp);
}

static double keygen_cloud(System* sys) {
    auto st = Clock::now();

    element_init_Zr(sys->skcs, sys->pairing);
    element_init_G1(sys->pkcs, sys->pairing);

    element_random(sys->skcs);
    element_pow_zn(sys->pkcs, sys->g, sys->skcs);

    auto ed = Clock::now();
    return elapsed_ms(st, ed);
}

static double keygen_clock(System* sys) {
    auto st = Clock::now();

    element_init_Zr(sys->skts, sys->pairing);
    element_init_G1(sys->pkts, sys->pairing);

    element_random(sys->skts);
    element_pow_zn(sys->pkts, sys->g, sys->skts);

    auto ed = Clock::now();
    return elapsed_ms(st, ed);
}

static double keygen_patient(System* sys) {
    auto st = Clock::now();

    element_init_Zr(sys->skp, sys->pairing);
    element_init_G1(sys->pkp, sys->pairing);

    element_random(sys->skp);
    element_pow_zn(sys->pkp, sys->g, sys->skp);

    auto ed = Clock::now();
    return elapsed_ms(st, ed);
}

static double keygen_doctor(System* sys, DoctorKey* dk) {
    auto st = Clock::now();

    pairing_ptr pairing = sys->pairing;

    element_init_G1(dk->D1, pairing);
    element_init_G1(dk->D2, pairing);
    element_init_Zr(dk->D3, pairing);
    element_init_Zr(dk->D4, pairing);
    element_init_G1(dk->D5, pairing);
    element_init_G1(dk->pkd, pairing);

    element_t xd, y, delta;
    element_t fx, exp1, exp2, denom, inv, tmp;

    element_init_Zr(xd, pairing);
    element_init_Zr(y, pairing);
    element_init_Zr(delta, pairing);
    element_init_Zr(fx, pairing);
    element_init_Zr(exp1, pairing);
    element_init_Zr(exp2, pairing);
    element_init_Zr(denom, pairing);
    element_init_Zr(inv, pairing);
    element_init_Zr(tmp, pairing);

    element_random(xd);
    element_random(y);
    element_random(delta);

    element_set(dk->D3, delta);
    element_set(dk->D4, y);

    // pkd=g^y
    element_pow_zn(dk->pkd, sys->g, y);

    // D1=g2^{f(xd)*(-x1)/(xd-x1)}
    poly_f(fx, sys, xd);
    element_sub(denom, xd, sys->x1);
    element_invert(inv, denom);
    element_neg(tmp, sys->x1);
    element_mul(exp1, fx, tmp);
    element_mul(exp1, exp1, inv);
    element_pow_zn(dk->D1, sys->g2, exp1);

    // D2=g2^{b1*(-xd)/(x1-xd)}
    element_sub(denom, sys->x1, xd);
    element_invert(inv, denom);
    element_neg(tmp, xd);
    element_mul(exp2, sys->b1, tmp);
    element_mul(exp2, exp2, inv);
    element_pow_zn(dk->D2, sys->g2, exp2);

    // D5=g^{(alpha-b1*b2)/(lambda+delta)}
    element_mul(tmp, sys->b1, sys->b2);
    element_sub(exp1, sys->alpha, tmp);
    element_add(denom, sys->lambda, delta);
    element_invert(inv, denom);
    element_mul(exp1, exp1, inv);
    element_pow_zn(dk->D5, sys->g, exp1);

    element_clear(xd);
    element_clear(y);
    element_clear(delta);
    element_clear(fx);
    element_clear(exp1);
    element_clear(exp2);
    element_clear(denom);
    element_clear(inv);
    element_clear(tmp);

    auto ed = Clock::now();
    return elapsed_ms(st, ed);
}

static double dvs_prove(System* sys, element_t sig, element_t pk_receiver) {
    auto st = Clock::now();

    // Gao 论文没有给 DVS 具体构造，这里用 pairing-based DVS 常见形式模拟生成
    element_t r, h, part1, part2;

    element_init_Zr(r, sys->pairing);
    element_init_Zr(h, sys->pairing);
    element_init_G1(part1, sys->pairing);
    element_init_G1(part2, sys->pairing);
    element_init_G1(sig, sys->pairing);

    element_random(r);
    element_random(h);

    element_pow_zn(part1, sys->g, r);
    element_pow_zn(part2, pk_receiver, h);
    element_mul(sig, part1, part2);

    element_clear(r);
    element_clear(h);
    element_clear(part1);
    element_clear(part2);

    auto ed = Clock::now();
    return elapsed_ms(st, ed);
}

static double dvs_verify(System* sys, element_t sig, element_t pk_receiver) {
    auto st = Clock::now();

    element_t h, tempG, gt1, gt2;

    element_init_Zr(h, sys->pairing);
    element_init_G1(tempG, sys->pairing);
    element_init_GT(gt1, sys->pairing);
    element_init_GT(gt2, sys->pairing);

    element_random(h);
    element_pow_zn(tempG, sys->g, h);

    pairing_apply(gt1, sig, sys->g, sys->pairing);
    pairing_apply(gt2, sys->pkts, pk_receiver, sys->pairing);

    element_clear(h);
    element_clear(tempG);
    element_clear(gt1);
    element_clear(gt2);

    auto ed = Clock::now();
    return elapsed_ms(st, ed);
}

static double ciphertext_generation(System* sys, Ciphertext* cw, int n) {
    auto st = Clock::now();

    pairing_ptr pairing = sys->pairing;

    cw->C1 = alloc_elements(n + 1);
    for (int i = 0; i <= n; i++) {
        element_init_G1(cw->C1[i], pairing);
    }

    element_init_G1(cw->C2, pairing);
    element_init_G1(cw->C3, pairing);

    element_t r, pi;
    element_init_Zr(r, pairing);
    element_init_Zr(pi, pairing);

    element_random(r);

    // C1,0=g2^r
    element_pow_zn(cw->C1[0], sys->g2, r);

    // C1,psi=g1^{r*pi_psi}
    for (int i = 1; i <= n; i++) {
        element_random(pi);
        element_mul(pi, pi, r);
        element_pow_zn(cw->C1[i], sys->g1, pi);
    }

    // C2=EK^r
    element_pow_zn(cw->C2, sys->EK, r);

    // C3=g^r
    element_pow_zn(cw->C3, sys->g, r);

    element_clear(r);
    element_clear(pi);

    auto ed = Clock::now();
    return elapsed_ms(st, ed);
}

static double trapdoor_generation(System* sys, DoctorKey* dk, Trapdoor* tq, int n, int m) {
    auto st = Clock::now();

    pairing_ptr pairing = sys->pairing;

    tq->T1 = alloc_elements(n + 1);
    for (int i = 0; i <= n; i++) {
        element_init_G1(tq->T1[i], pairing);
    }

    element_init_G1(tq->T2, pairing);
    element_init_G1(tq->T3, pairing);
    element_init_G1(tq->T4, pairing);
    element_init_Zr(tq->T5, pairing);

    element_t s, tau, exp, sum, hpow, temp;

    element_init_Zr(s, pairing);
    element_init_Zr(tau, pairing);
    element_init_Zr(exp, pairing);
    element_init_Zr(sum, pairing);
    element_init_Zr(hpow, pairing);
    element_init_G1(temp, pairing);

    element_random(s);
    element_random(tau);

    // T1,0=g1^s
    element_pow_zn(tq->T1[0], sys->g1, s);

    // T1,psi=g2^{s * sum H(w)^psi} * pkcs^tau
    for (int psi = 1; psi <= n; psi++) {
        element_set0(sum);

        for (int q = 0; q < m; q++) {
            element_random(hpow);
            element_add(sum, sum, hpow);
        }

        element_mul(exp, s, sum);
        element_pow_zn(tq->T1[psi], sys->g2, exp);

        element_pow_zn(temp, sys->pkcs, tau);
        element_mul(tq->T1[psi], tq->T1[psi], temp);
    }

    // T2=D2^s
    element_pow_zn(tq->T2, dk->D2, s);

    // T3=D1^s
    element_pow_zn(tq->T3, dk->D1, s);

    // T4=g^tau
    element_pow_zn(tq->T4, sys->g, tau);

    // T5=D3
    element_set(tq->T5, dk->D3);

    element_clear(s);
    element_clear(tau);
    element_clear(exp);
    element_clear(sum);
    element_clear(hpow);
    element_clear(temp);

    auto ed = Clock::now();
    return elapsed_ms(st, ed);
}

static double ciphertext_match(System* sys, Ciphertext* cw, Trapdoor* tq, int n) {
    auto st = Clock::now();

    pairing_ptr pairing = sys->pairing;

    element_t left, right, tempGT, tempG, inv;

    element_init_GT(left, pairing);
    element_init_GT(right, pairing);
    element_init_GT(tempGT, pairing);
    element_init_G1(tempG, pairing);
    element_init_G1(inv, pairing);

    // e(C1,0,T1,0)
    pairing_apply(left, cw->C1[0], tq->T1[0], pairing);

    // prod e(C1,psi, T1,psi / T4^skcs)
    for (int i = 1; i <= n; i++) {
        element_pow_zn(tempG, tq->T4, sys->skcs);
        element_invert(inv, tempG);
        element_mul(tempG, tq->T1[i], inv);

        pairing_apply(tempGT, cw->C1[i], tempG, pairing);
        element_mul(left, left, tempGT);
    }

    // e(C2,T2) * e(C3,T3)
    pairing_apply(right, cw->C2, tq->T2, pairing);
    pairing_apply(tempGT, cw->C3, tq->T3, pairing);
    element_mul(right, right, tempGT);

    element_clear(left);
    element_clear(right);
    element_clear(tempGT);
    element_clear(tempG);
    element_clear(inv);

    auto ed = Clock::now();
    return elapsed_ms(st, ed);
}

static void clear_ciphertext(Ciphertext* cw, int n) {
    for (int i = 0; i <= n; i++) {
        element_clear(cw->C1[i]);
    }
    std::free(cw->C1);
    element_clear(cw->C2);
    element_clear(cw->C3);
}

static void clear_trapdoor(Trapdoor* tq, int n) {
    for (int i = 0; i <= n; i++) {
        element_clear(tq->T1[i]);
    }
    std::free(tq->T1);
    element_clear(tq->T2);
    element_clear(tq->T3);
    element_clear(tq->T4);
    element_clear(tq->T5);
}

static void clear_doctor_key(DoctorKey* dk) {
    element_clear(dk->D1);
    element_clear(dk->D2);
    element_clear(dk->D3);
    element_clear(dk->D4);
    element_clear(dk->D5);
    element_clear(dk->pkd);
}

static void clear_system(System* sys) {
    element_clear(sys->g);
    element_clear(sys->g1);
    element_clear(sys->g2);
    element_clear(sys->Y);
    element_clear(sys->EK);

    element_clear(sys->b1);
    element_clear(sys->b2);
    element_clear(sys->a1);
    element_clear(sys->x1);
    element_clear(sys->alpha);
    element_clear(sys->lambda);

    element_clear(sys->skcs);
    element_clear(sys->pkcs);
    element_clear(sys->skts);
    element_clear(sys->pkts);
    element_clear(sys->skp);
    element_clear(sys->pkp);
}

int main(int argc, char** argv) {
    Params p = { 10, 2 };

    if (argc > 1) p.n = std::atoi(argv[1]);
    if (argc > 2) p.m = std::atoi(argv[2]);

    pairing_t pairing;
    pairing_init_set_buf(pairing, DEFAULT_PAIRING_PARAM, std::strlen(DEFAULT_PAIRING_PARAM));

    std::printf("blocks,dec_keygen_ms,sign_keygen_ms,encryption_ms,signature_ms,decryption_ms,verification_ms\n");

    for (int blocks = 50; blocks <= 1000; blocks += 50) {
        System sys;
        DoctorKey dk;

        init_system(&sys, pairing);

        // 辅助密钥生成：用于算法运行，但不计入六项输出
        keygen_cloud(&sys);
        keygen_patient(&sys);

        // 六项中的签名密钥生成：KeyGents
        double t_sign_keygen = keygen_clock(&sys);

        // 六项中的解密密钥生成：KeyGend
        double t_dec_keygen = keygen_doctor(&sys, &dk);

        double t_encryption = 0.0;
        double t_signature = 0.0;
        double t_decryption = 0.0;
        double t_verification = 0.0;

        for (int b = 0; b < blocks; b++) {
            Ciphertext cw;
            Trapdoor tq;
            element_t SPT, SDT;

            // Patient Time + Doctor Time 里的 DVS.Prove
            t_signature += dvs_prove(&sys, SPT, sys.pkp);
            t_signature += dvs_prove(&sys, SDT, dk.pkd);

            // EHR index ciphertext generation
            t_encryption += ciphertext_generation(&sys, &cw, p.n);

            // Trapdoor generation + ciphertext match 归入解密/检索时间
            t_decryption += trapdoor_generation(&sys, &dk, &tq, p.n, p.m);
            t_decryption += ciphertext_match(&sys, &cw, &tq, p.n);

            // Time Match 中的 DVS.Verify
            t_verification += dvs_verify(&sys, SPT, sys.pkp);
            t_verification += dvs_verify(&sys, SDT, dk.pkd);

            element_clear(SPT);
            element_clear(SDT);
            clear_trapdoor(&tq, p.n);
            clear_ciphertext(&cw, p.n);
        }

        std::printf("%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            blocks,
            t_dec_keygen,
            t_sign_keygen,
            t_encryption,
            t_signature,
            t_decryption,
            t_verification);

        clear_doctor_key(&dk);
        clear_system(&sys);
    }

    pairing_clear(pairing);
    return 0;
}