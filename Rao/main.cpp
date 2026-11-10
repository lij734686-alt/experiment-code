#define _CRT_SECURE_NO_WARNINGS
#include <pbc.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using Clock = std::chrono::high_resolution_clock;

// 为了让时间短一点，把属性和关键词规模设小
static const int SIGN_ATTRS = 5;
static const int DEC_ATTRS = 5;
static const int KW_ROWS = 2;
static const int CLAUSE_ATTRS = 3;

static const char* DEFAULT_PARAM =
"type a\n"
"q 8780710799663312522437781984754049815806883199414208211028653399266475630880222957078625179422662221423155858769582317459277713367317481324925129998224791\n"
"h 120160122648911460793888213667405342048029544012513118229196151310472072893597045311\n"
"r 730750818665451621361119245571504901405976559617\n"
"exp2 159\n"
"exp1 107\n"
"sign1 1\n"
"sign0 1\n";

static volatile int g_sink = 0;

double ms_since(Clock::time_point st, Clock::time_point ed) {
    return std::chrono::duration<double, std::milli>(ed - st).count();
}

void H_G1(element_t out, const char* tag, int i, int j) {
    char buf[128];
    std::snprintf(buf, sizeof(buf), "%s_%d_%d", tag, i, j);
    element_from_hash(out, buf, (int)std::strlen(buf));
}

void H_Zr(element_t out, const char* tag, int i, int j) {
    char buf[128];
    std::snprintf(buf, sizeof(buf), "%s_%d_%d", tag, i, j);
    element_from_hash(out, buf, (int)std::strlen(buf));
}

struct System {
    pairing_t pairing;
    element_t alpha, z1, z2, z3, beta;
    element_t g, h, g1, g2, g3, g4, g5;
    element_t h1, h2, h3;
    element_t gT, Y;
};

struct SignKey {
    element_t S01, S02;
    element_t Sx[SIGN_ATTRS];
};

struct DecKey {
    element_t D01, D02, D03;
    element_t Dy[DEC_ATTRS];
};

struct Ciphertext {
    element_t theta, gamma, delta, eta;
    element_t K, K1, K2;
    element_t KW[KW_ROWS], LW[KW_ROWS];
    element_t E, Ei1, Ei2, E0, Eprime;
    element_t sigma, sigma_i[SIGN_ATTRS];
    element_t tage, tags;
};

struct Token {
    element_t lambda, tau;
    element_t T1, T2, tk;
    element_t Ti1[KW_ROWS], Ti2[KW_ROWS];
    element_t D01p, D02p, D02pp, D03p;
    element_t Dyp[DEC_ATTRS];
};

struct MatchCT {
    element_t X1;
};

struct TransformCT {
    element_t X1, X2;
};

void init_system(System& s) {
    pairing_init_set_buf(s.pairing, DEFAULT_PARAM, std::strlen(DEFAULT_PARAM));

    element_init_Zr(s.alpha, s.pairing);
    element_init_Zr(s.z1, s.pairing);
    element_init_Zr(s.z2, s.pairing);
    element_init_Zr(s.z3, s.pairing);
    element_init_Zr(s.beta, s.pairing);

    element_init_G1(s.g, s.pairing);
    element_init_G1(s.h, s.pairing);
    element_init_G1(s.g1, s.pairing);
    element_init_G1(s.g2, s.pairing);
    element_init_G1(s.g3, s.pairing);
    element_init_G1(s.g4, s.pairing);
    element_init_G1(s.g5, s.pairing);
    element_init_G1(s.h1, s.pairing);
    element_init_G1(s.h2, s.pairing);
    element_init_G1(s.h3, s.pairing);
    element_init_G1(s.Y, s.pairing);
    element_init_GT(s.gT, s.pairing);

    element_random(s.alpha);
    element_random(s.z1);
    element_random(s.z2);
    element_random(s.z3);
    element_random(s.beta);

    element_random(s.g);
    element_random(s.h);
    element_random(s.g1);
    element_random(s.g2);
    element_random(s.g3);
    element_random(s.g4);
    element_random(s.g5);

    element_t tmp;
    element_init_GT(tmp, s.pairing);
    element_pairing(tmp, s.g, s.g);
    element_pow_zn(s.gT, tmp, s.alpha);

    element_pow_zn(s.h1, s.g, s.z1);
    element_pow_zn(s.h2, s.g, s.z2);
    element_pow_zn(s.h3, s.g, s.z3);
    element_pow_zn(s.Y, s.g, s.beta);

    element_clear(tmp);
}

void clear_system(System& s) {
    element_clear(s.alpha);
    element_clear(s.z1);
    element_clear(s.z2);
    element_clear(s.z3);
    element_clear(s.beta);
    element_clear(s.g);
    element_clear(s.h);
    element_clear(s.g1);
    element_clear(s.g2);
    element_clear(s.g3);
    element_clear(s.g4);
    element_clear(s.g5);
    element_clear(s.h1);
    element_clear(s.h2);
    element_clear(s.h3);
    element_clear(s.Y);
    element_clear(s.gT);
    pairing_clear(s.pairing);
}

void init_sign_key(SignKey& sk, System& s) {
    element_init_G1(sk.S01, s.pairing);
    element_init_G1(sk.S02, s.pairing);
    for (int i = 0; i < SIGN_ATTRS; i++) element_init_G1(sk.Sx[i], s.pairing);
}

void clear_sign_key(SignKey& sk) {
    element_clear(sk.S01);
    element_clear(sk.S02);
    for (int i = 0; i < SIGN_ATTRS; i++) element_clear(sk.Sx[i]);
}

void init_dec_key(DecKey& dk, System& s) {
    element_init_G1(dk.D01, s.pairing);
    element_init_G1(dk.D02, s.pairing);
    element_init_G1(dk.D03, s.pairing);
    for (int i = 0; i < DEC_ATTRS; i++) element_init_G1(dk.Dy[i], s.pairing);
}

void clear_dec_key(DecKey& dk) {
    element_clear(dk.D01);
    element_clear(dk.D02);
    element_clear(dk.D03);
    for (int i = 0; i < DEC_ATTRS; i++) element_clear(dk.Dy[i]);
}

void init_ct(Ciphertext& ct, System& s) {
    element_init_Zr(ct.theta, s.pairing);
    element_init_Zr(ct.gamma, s.pairing);
    element_init_Zr(ct.delta, s.pairing);
    element_init_Zr(ct.eta, s.pairing);

    element_init_G1(ct.K, s.pairing);
    element_init_G1(ct.K1, s.pairing);
    element_init_G1(ct.K2, s.pairing);
    element_init_G1(ct.E, s.pairing);
    element_init_G1(ct.Ei1, s.pairing);
    element_init_G1(ct.Ei2, s.pairing);
    element_init_G1(ct.E0, s.pairing);
    element_init_G1(ct.Eprime, s.pairing);
    element_init_G1(ct.sigma, s.pairing);

    for (int i = 0; i < KW_ROWS; i++) {
        element_init_G1(ct.KW[i], s.pairing);
        element_init_G1(ct.LW[i], s.pairing);
    }

    for (int i = 0; i < SIGN_ATTRS; i++) {
        element_init_G1(ct.sigma_i[i], s.pairing);
    }

    element_init_Zr(ct.tage, s.pairing);
    element_init_Zr(ct.tags, s.pairing);
}

void clear_ct(Ciphertext& ct) {
    element_clear(ct.theta);
    element_clear(ct.gamma);
    element_clear(ct.delta);
    element_clear(ct.eta);
    element_clear(ct.K);
    element_clear(ct.K1);
    element_clear(ct.K2);
    element_clear(ct.E);
    element_clear(ct.Ei1);
    element_clear(ct.Ei2);
    element_clear(ct.E0);
    element_clear(ct.Eprime);
    element_clear(ct.sigma);

    for (int i = 0; i < KW_ROWS; i++) {
        element_clear(ct.KW[i]);
        element_clear(ct.LW[i]);
    }

    for (int i = 0; i < SIGN_ATTRS; i++) {
        element_clear(ct.sigma_i[i]);
    }

    element_clear(ct.tage);
    element_clear(ct.tags);
}

void init_token(Token& tk, System& s) {
    element_init_Zr(tk.lambda, s.pairing);
    element_init_Zr(tk.tau, s.pairing);
    element_init_G1(tk.T1, s.pairing);
    element_init_G1(tk.T2, s.pairing);
    element_init_GT(tk.tk, s.pairing);

    for (int i = 0; i < KW_ROWS; i++) {
        element_init_G1(tk.Ti1[i], s.pairing);
        element_init_G1(tk.Ti2[i], s.pairing);
    }

    element_init_G1(tk.D01p, s.pairing);
    element_init_G1(tk.D02p, s.pairing);
    element_init_G1(tk.D02pp, s.pairing);
    element_init_G1(tk.D03p, s.pairing);

    for (int i = 0; i < DEC_ATTRS; i++) {
        element_init_G1(tk.Dyp[i], s.pairing);
    }
}

void clear_token(Token& tk) {
    element_clear(tk.lambda);
    element_clear(tk.tau);
    element_clear(tk.T1);
    element_clear(tk.T2);
    element_clear(tk.tk);

    for (int i = 0; i < KW_ROWS; i++) {
        element_clear(tk.Ti1[i]);
        element_clear(tk.Ti2[i]);
    }

    element_clear(tk.D01p);
    element_clear(tk.D02p);
    element_clear(tk.D02pp);
    element_clear(tk.D03p);

    for (int i = 0; i < DEC_ATTRS; i++) {
        element_clear(tk.Dyp[i]);
    }
}

void init_match(MatchCT& m, System& s) {
    element_init_GT(m.X1, s.pairing);
}

void clear_match(MatchCT& m) {
    element_clear(m.X1);
}

void init_tr(TransformCT& tr, System& s) {
    element_init_GT(tr.X1, s.pairing);
    element_init_GT(tr.X2, s.pairing);
}

void clear_tr(TransformCT& tr) {
    element_clear(tr.X1);
    element_clear(tr.X2);
}

// 签名密钥生成：S01 = g^alpha h^rs, S02 = g^rs, Sx = H1(x)^rs
void s_keygen(System& s, SignKey& sk) {
    element_t rs, tmp, hx;
    element_init_Zr(rs, s.pairing);
    element_init_G1(tmp, s.pairing);
    element_init_G1(hx, s.pairing);

    element_random(rs);

    element_pow_zn(sk.S01, s.g, s.alpha);
    element_pow_zn(tmp, s.h, rs);
    element_mul(sk.S01, sk.S01, tmp);

    element_pow_zn(sk.S02, s.g, rs);

    for (int i = 0; i < SIGN_ATTRS; i++) {
        H_G1(hx, "H1_SIGN_ATTR", i, 0);
        element_pow_zn(sk.Sx[i], hx, rs);
    }

    element_clear(rs);
    element_clear(tmp);
    element_clear(hx);
}

// 解密密钥生成
void d_keygen(System& s, DecKey& dk) {
    element_t re, t1, t2, exp, tmp, hy;
    element_init_Zr(re, s.pairing);
    element_init_Zr(t1, s.pairing);
    element_init_Zr(t2, s.pairing);
    element_init_Zr(exp, s.pairing);
    element_init_G1(tmp, s.pairing);
    element_init_G1(hy, s.pairing);

    element_random(re);

    element_mul(t1, s.z1, s.z3);
    element_mul(t1, t1, re);
    element_add(exp, s.alpha, t1);
    element_pow_zn(dk.D01, s.g, exp);
    element_pow_zn(tmp, s.h, re);
    element_mul(dk.D01, dk.D01, tmp);

    element_pow_zn(dk.D02, s.g, re);

    element_mul(t1, s.z1, s.z3);
    element_sub(t1, t1, re);
    element_invert(t2, s.z2);
    element_mul(exp, t1, t2);
    element_pow_zn(dk.D03, s.g, exp);

    for (int i = 0; i < DEC_ATTRS; i++) {
        H_G1(hy, "H2_DEC_ATTR", i, 0);
        element_pow_zn(dk.Dy[i], hy, re);
    }

    element_clear(re);
    element_clear(t1);
    element_clear(t2);
    element_clear(exp);
    element_clear(tmp);
    element_clear(hy);
}

// Signcrypt 的加密部分
void signcrypt_encrypt_part(System& s, Ciphertext& ct, int block_id) {
    element_t o1, o2, h9, inv_h9, exp, tmpG, tmpG2, h8, wj, piw;
    element_t pairTmp;

    element_init_Zr(o1, s.pairing);
    element_init_Zr(o2, s.pairing);
    element_init_Zr(h9, s.pairing);
    element_init_Zr(inv_h9, s.pairing);
    element_init_Zr(exp, s.pairing);
    element_init_Zr(wj, s.pairing);
    element_init_Zr(piw, s.pairing);
    element_init_G1(tmpG, s.pairing);
    element_init_G1(tmpG2, s.pairing);
    element_init_G1(h8, s.pairing);
    element_init_GT(pairTmp, s.pairing);

    element_random(ct.theta);
    H_Zr(ct.gamma, "H3_KDF_KEY", block_id, 0);
    H_Zr(ct.delta, "H7_DELTA", block_id, 0);
    element_random(ct.eta);

    element_random(o1);
    element_random(o2);
    element_pow_zn(ct.K1, s.g, o1);
    element_pow_zn(ct.K2, s.g, o2);
    element_pairing(pairTmp, ct.K1, s.Y);
    element_pow_zn(pairTmp, pairTmp, o2);

    H_G1(h8, "H8_ENC", block_id, 0);
    element_pow_zn(tmpG, s.g, ct.theta);
    element_mul(ct.E, h8, tmpG);

    element_random(exp);
    element_pow_zn(ct.Ei1, s.g, exp);

    element_pow_zn(ct.Ei2, s.h, exp);
    for (int i = 0; i < CLAUSE_ATTRS; i++) {
        H_G1(tmpG, "H2_DEC_ATTR", i, 0);
        element_pow_zn(tmpG2, tmpG, exp);
        element_mul(ct.Ei2, ct.Ei2, tmpG2);
    }

    element_pow_zn(ct.Eprime, s.g, ct.theta);

    H_Zr(h9, "H9_HT", block_id, 0);
    element_invert(inv_h9, h9);
    element_mul(exp, ct.theta, inv_h9);
    element_pow_zn(ct.K, s.h2, exp);

    for (int j = 0; j < KW_ROWS; j++) {
        element_random(piw);

        // 与 token_gen 保持一致，模拟同一组关键词
        H_Zr(wj, "KEYWORD_VALUE", j, 0);

        element_add(exp, ct.theta, piw);
        element_pow_zn(ct.KW[j], s.h1, exp);

        element_mul(exp, wj, piw);
        element_pow_zn(tmpG, s.h2, exp);
        element_mul(ct.KW[j], ct.KW[j], tmpG);

        element_pow_zn(ct.LW[j], s.h3, piw);
    }

    element_pow_zn(tmpG, s.g3, ct.delta);
    element_pow_zn(tmpG2, s.g4, ct.eta);
    element_mul(tmpG, tmpG, tmpG2);
    element_mul(tmpG, tmpG, s.g5);
    element_pow_zn(ct.E0, tmpG, ct.theta);

    H_Zr(ct.tage, "H4_TAGE", block_id, 0);
    H_Zr(ct.tags, "H5_TAGS", block_id, 0);

    element_clear(o1);
    element_clear(o2);
    element_clear(h9);
    element_clear(inv_h9);
    element_clear(exp);
    element_clear(wj);
    element_clear(piw);
    element_clear(tmpG);
    element_clear(tmpG2);
    element_clear(h8);
    element_clear(pairTmp);
}

// Signcrypt 的签名部分
void signcrypt_signature_part(System& s, SignKey& sk, Ciphertext& ct, int block_id) {
    element_t rs2, inv_gamma, o2, o_sig, sum_b, exp;
    element_t tmpG, tmpG2, hAttr, S01p, S02p, base, inner;
    element_t b[SIGN_ATTRS];

    element_init_Zr(rs2, s.pairing);
    element_init_Zr(inv_gamma, s.pairing);
    element_init_Zr(o2, s.pairing);
    element_init_Zr(o_sig, s.pairing);
    element_init_Zr(sum_b, s.pairing);
    element_init_Zr(exp, s.pairing);
    element_init_G1(tmpG, s.pairing);
    element_init_G1(tmpG2, s.pairing);
    element_init_G1(hAttr, s.pairing);
    element_init_G1(S01p, s.pairing);
    element_init_G1(S02p, s.pairing);
    element_init_G1(base, s.pairing);
    element_init_G1(inner, s.pairing);

    for (int i = 0; i < SIGN_ATTRS; i++) element_init_Zr(b[i], s.pairing);

    element_random(rs2);
    element_random(o_sig);
    H_Zr(o2, "H6_SIGN_HASH", block_id, 0);
    element_invert(inv_gamma, ct.gamma);

    element_pow_zn(tmpG, s.h, rs2);
    element_mul(S01p, sk.S01, tmpG);

    element_pow_zn(tmpG, s.g, rs2);
    element_mul(S02p, sk.S02, tmpG);

    element_set0(sum_b);
    for (int i = 0; i < SIGN_ATTRS - 1; i++) {
        element_random(b[i]);
        element_add(sum_b, sum_b, b[i]);
    }
    element_neg(b[SIGN_ATTRS - 1], sum_b);

    element_pow_zn(tmpG, s.g1, o2);
    element_mul(base, tmpG, s.g2);
    element_pow_zn(base, base, ct.theta);

    element_mul(inner, S01p, base);
    element_pow_zn(ct.sigma, inner, inv_gamma);

    for (int i = 0; i < SIGN_ATTRS; i++) {
        H_G1(hAttr, "H1_SIGN_ATTR", i, 0);

        element_pow_zn(tmpG, sk.Sx[i], inv_gamma);

        element_mul(exp, o_sig, b[i]);
        element_pow_zn(tmpG2, hAttr, exp);

        element_mul(tmpG, tmpG, tmpG2);
        element_mul(ct.sigma, ct.sigma, tmpG);

        element_pow_zn(ct.sigma_i[i], S02p, inv_gamma);
        element_pow_zn(tmpG, s.g, exp);
        element_mul(ct.sigma_i[i], ct.sigma_i[i], tmpG);
    }

    for (int i = 0; i < SIGN_ATTRS; i++) element_clear(b[i]);

    element_clear(rs2);
    element_clear(inv_gamma);
    element_clear(o2);
    element_clear(o_sig);
    element_clear(sum_b);
    element_clear(exp);
    element_clear(tmpG);
    element_clear(tmpG2);
    element_clear(hAttr);
    element_clear(S01p);
    element_clear(S02p);
    element_clear(base);
    element_clear(inner);
}

// TokenGen：优化后只生成一次
void token_gen(System& s, DecKey& dk, Token& tk) {
    element_t f1, f2, inv_lambda, inv_tau, mu, wj;
    element_t tmpG, h8, pairTmp;

    element_init_Zr(f1, s.pairing);
    element_init_Zr(f2, s.pairing);
    element_init_Zr(inv_lambda, s.pairing);
    element_init_Zr(inv_tau, s.pairing);
    element_init_Zr(mu, s.pairing);
    element_init_Zr(wj, s.pairing);
    element_init_G1(tmpG, s.pairing);
    element_init_G1(h8, s.pairing);
    element_init_GT(pairTmp, s.pairing);

    element_random(f1);
    element_random(f2);
    element_random(tk.lambda);
    element_random(tk.tau);

    element_invert(inv_lambda, tk.lambda);
    element_invert(inv_tau, tk.tau);

    element_pow_zn(tk.T1, s.g, f1);
    element_pow_zn(tk.T2, s.g, f2);
    element_pairing(pairTmp, tk.T1, s.Y);
    element_pow_zn(tk.tk, pairTmp, f2);

    for (int i = 0; i < KW_ROWS; i++) {
        element_random(mu);

        H_Zr(wj, "KEYWORD_VALUE", i, 0);

        element_pow_zn(tmpG, s.h2, wj);
        element_mul(tmpG, s.h1, tmpG);
        element_pow_zn(tk.Ti1[i], tmpG, mu);

        H_G1(h8, "H8_TOKEN_ROW", i, 0);
        element_pow_zn(tmpG, s.h3, mu);
        element_mul(tk.Ti2[i], h8, tmpG);
    }

    element_pow_zn(tk.D01p, dk.D01, inv_tau);
    element_pow_zn(tk.D02p, dk.D02, inv_lambda);
    element_pow_zn(tk.D02pp, dk.D02, inv_tau);
    element_pow_zn(tk.D03p, dk.D03, inv_lambda);

    for (int i = 0; i < DEC_ATTRS; i++) {
        element_pow_zn(tk.Dyp[i], dk.Dy[i], inv_tau);
    }

    element_clear(f1);
    element_clear(f2);
    element_clear(inv_lambda);
    element_clear(inv_tau);
    element_clear(mu);
    element_clear(wj);
    element_clear(tmpG);
    element_clear(h8);
    element_clear(pairTmp);
}

void search_phase(System& s, Ciphertext& ct, Token& tk, MatchCT& m, int block_id) {
    element_t p1, p2, tmpGT, invGT;
    element_t tmpG, h9;

    element_init_GT(p1, s.pairing);
    element_init_GT(p2, s.pairing);
    element_init_GT(tmpGT, s.pairing);
    element_init_GT(invGT, s.pairing);
    element_init_G1(tmpG, s.pairing);
    element_init_Zr(h9, s.pairing);

    element_set1(m.X1);

    for (int i = 0; i < KW_ROWS; i++) {
        element_pairing(p1, ct.KW[i], tk.Ti2[i]);
        element_pairing(p2, ct.LW[i], tk.Ti1[i]);
        element_invert(invGT, p2);
        element_mul(tmpGT, p1, invGT);
        element_mul(m.X1, m.X1, tmpGT);
    }

    H_Zr(h9, "H9_HT", block_id, 0);
    element_pow_zn(tmpG, ct.K, h9);
    element_pairing(p1, tk.D02p, ct.Eprime);
    element_pairing(p2, tk.D03p, tmpG);
    element_mul(tmpGT, p1, p2);

    g_sink ^= element_is1(tmpGT);

    element_clear(p1);
    element_clear(p2);
    element_clear(tmpGT);
    element_clear(invGT);
    element_clear(tmpG);
    element_clear(h9);
}

void transform_phase(System& s, Ciphertext& ct, Token& tk, MatchCT& m, TransformCT& tr) {
    element_t prodDy, p1, p2, p3, invGT;

    element_init_G1(prodDy, s.pairing);
    element_init_GT(p1, s.pairing);
    element_init_GT(p2, s.pairing);
    element_init_GT(p3, s.pairing);
    element_init_GT(invGT, s.pairing);

    element_set1(prodDy);
    for (int i = 0; i < CLAUSE_ATTRS; i++) {
        element_mul(prodDy, prodDy, tk.Dyp[i]);
    }

    element_pairing(p1, ct.Ei2, tk.D02pp);
    element_invert(invGT, p1);

    element_pairing(p2, ct.Ei1, prodDy);
    element_pairing(p3, ct.Eprime, tk.D01p);

    element_mul(tr.X2, invGT, p2);
    element_mul(tr.X2, tr.X2, p3);
    element_set(tr.X1, m.X1);

    element_clear(prodDy);
    element_clear(p1);
    element_clear(p2);
    element_clear(p3);
    element_clear(invGT);
}

void unsigncrypt_recover(System& s, Token& tk, TransformCT& tr, element_t J) {
    element_t neg_lambda, a, b;

    element_init_Zr(neg_lambda, s.pairing);
    element_init_GT(a, s.pairing);
    element_init_GT(b, s.pairing);

    element_neg(neg_lambda, tk.lambda);
    element_pow_zn(a, tr.X1, neg_lambda);
    element_pow_zn(b, tr.X2, tk.tau);
    element_mul(J, a, b);

    element_clear(neg_lambda);
    element_clear(a);
    element_clear(b);
}

void verify_phase(System& s, Ciphertext& ct, TransformCT& tr, element_t J, int block_id) {
    element_t gamma_check, inv_gamma, a, b, c;

    element_init_Zr(gamma_check, s.pairing);
    element_init_Zr(inv_gamma, s.pairing);
    element_init_GT(a, s.pairing);
    element_init_GT(b, s.pairing);
    element_init_GT(c, s.pairing);

    H_Zr(gamma_check, "H3_VERIFY", block_id, 0);
    element_invert(inv_gamma, gamma_check);

    element_pow_zn(a, J, gamma_check);
    H_Zr(ct.tage, "H4_VERIFY_TAGE", block_id, 0);

    element_pow_zn(b, s.gT, inv_gamma);
    H_Zr(ct.tags, "H5_VERIFY_TAGS", block_id, 0);

    element_mul(c, a, b);
    g_sink ^= element_is1(c);

    element_clear(gamma_check);
    element_clear(inv_gamma);
    element_clear(a);
    element_clear(b);
    element_clear(c);
}

int main() {
    std::printf("blocks,dec_keygen_ms,sign_keygen_ms,encryption_ms,signature_ms,decryption_ms,verification_ms\n");

    for (int blocks = 50; blocks <= 1000; blocks += 50) {
        System sys;
        init_system(sys);

        SignKey sk;
        DecKey dk;
        init_sign_key(sk, sys);
        init_dec_key(dk, sys);

        auto st = Clock::now();
        d_keygen(sys, dk);
        auto ed = Clock::now();
        double dec_keygen_ms = ms_since(st, ed);

        st = Clock::now();
        s_keygen(sys, sk);
        ed = Clock::now();
        double sign_keygen_ms = ms_since(st, ed);

        double encryption_ms = 0.0;
        double signature_ms = 0.0;
        double decryption_ms = 0.0;
        double verification_ms = 0.0;

        // 固定 1 个用户和 1 次查询，TokenGen 只执行一次
        Token tk;
        init_token(tk, sys);

        st = Clock::now();
        token_gen(sys, dk, tk);
        ed = Clock::now();
        decryption_ms += ms_since(st, ed);

        for (int b = 0; b < blocks; b++) {
            Ciphertext ct;
            MatchCT mt;
            TransformCT tr;

            init_ct(ct, sys);
            init_match(mt, sys);
            init_tr(tr, sys);

            element_t J;
            element_init_GT(J, sys.pairing);

            st = Clock::now();
            signcrypt_encrypt_part(sys, ct, b);
            ed = Clock::now();
            encryption_ms += ms_since(st, ed);

            st = Clock::now();
            signcrypt_signature_part(sys, sk, ct, b);
            ed = Clock::now();
            signature_ms += ms_since(st, ed);

            st = Clock::now();
            search_phase(sys, ct, tk, mt, b);
            transform_phase(sys, ct, tk, mt, tr);
            unsigncrypt_recover(sys, tk, tr, J);
            ed = Clock::now();
            decryption_ms += ms_since(st, ed);

            st = Clock::now();
            verify_phase(sys, ct, tr, J, b);
            ed = Clock::now();
            verification_ms += ms_since(st, ed);

            element_clear(J);
            clear_tr(tr);
            clear_match(mt);
            clear_ct(ct);
        }

        std::printf("%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            blocks,
            dec_keygen_ms,
            sign_keygen_ms,
            encryption_ms,
            signature_ms,
            decryption_ms,
            verification_ms
        );

        clear_token(tk);
        clear_dec_key(dk);
        clear_sign_key(sk);
        clear_system(sys);
    }

    if (g_sink == 123456789) {
        std::printf("ignore: %d\n", g_sink);
    }

    return 0;
}