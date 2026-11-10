#define _CRT_SECURE_NO_WARNINGS
#include <pbc.h>
#include <chrono>
#include <cstdio>
#include <cstring>

using Clock = std::chrono::high_resolution_clock;

static const int ATTRS = 5;       // DU 属性数量
static const int S_TAG = 5;       // 标签中向量长度 s
static volatile int g_sink = 0;

static const char* DEFAULT_PARAM =
"type a\n"
"q 8780710799663312522437781984754049815806883199414208211028653399266475630880222957078625179422662221423155858769582317459277713367317481324925129998224791\n"
"h 120160122648911460793888213667405342048029544012513118229196151310472072893597045311\n"
"r 730750818665451621361119245571504901405976559617\n"
"exp2 159\n"
"exp1 107\n"
"sign1 1\n"
"sign0 1\n";

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

    element_t g;
    element_t beta;
    element_t alpha;
    element_t s;

    element_t g_beta;
    element_t MPK;
    element_t egg_alpha;
};

struct OwnerKey {
    element_t Kdos;
    element_t alphaDO;
    element_t Kdp;
};

struct AuditorKey {
    element_t Kdas;
    element_t alphaDA;
    element_t Ktp;
};

struct UserDecKey {
    element_t D1;
    element_t Dm[ATTRS];
    element_t Dm_p[ATTRS];
};

struct CipherBlock {
    element_t T1;
    element_t Ci;
    element_t c;
    element_t C1;
    element_t Cy[ATTRS];
    element_t Cy_p[ATTRS];
};

struct TagBlock {
    element_t sigma;
    element_t cij[S_TAG];
};

struct Proof {
    element_t sigma_agg;
    element_t mu[S_TAG];
};

void init_system(System& sys) {
    pairing_init_set_buf(sys.pairing, DEFAULT_PARAM, std::strlen(DEFAULT_PARAM));

    element_init_G1(sys.g, sys.pairing);
    element_init_Zr(sys.beta, sys.pairing);
    element_init_Zr(sys.alpha, sys.pairing);
    element_init_Zr(sys.s, sys.pairing);
    element_init_G1(sys.g_beta, sys.pairing);
    element_init_G1(sys.MPK, sys.pairing);
    element_init_GT(sys.egg_alpha, sys.pairing);

    element_random(sys.g);
    element_random(sys.beta);
    element_random(sys.alpha);
    element_random(sys.s);

    element_pow_zn(sys.g_beta, sys.g, sys.beta);
    element_pow_zn(sys.MPK, sys.g, sys.s);

    element_t egg;
    element_init_GT(egg, sys.pairing);
    element_pairing(egg, sys.g, sys.g);
    element_pow_zn(sys.egg_alpha, egg, sys.alpha);
    element_clear(egg);
}

void clear_system(System& sys) {
    element_clear(sys.g);
    element_clear(sys.beta);
    element_clear(sys.alpha);
    element_clear(sys.s);
    element_clear(sys.g_beta);
    element_clear(sys.MPK);
    element_clear(sys.egg_alpha);
    pairing_clear(sys.pairing);
}

void init_owner_key(OwnerKey& ok, System& sys) {
    element_init_G1(ok.Kdos, sys.pairing);
    element_init_Zr(ok.alphaDO, sys.pairing);
    element_init_G1(ok.Kdp, sys.pairing);
}

void clear_owner_key(OwnerKey& ok) {
    element_clear(ok.Kdos);
    element_clear(ok.alphaDO);
    element_clear(ok.Kdp);
}

void init_auditor_key(AuditorKey& ak, System& sys) {
    element_init_G1(ak.Kdas, sys.pairing);
    element_init_Zr(ak.alphaDA, sys.pairing);
    element_init_G1(ak.Ktp, sys.pairing);
}

void clear_auditor_key(AuditorKey& ak) {
    element_clear(ak.Kdas);
    element_clear(ak.alphaDA);
    element_clear(ak.Ktp);
}

void init_user_dec_key(UserDecKey& dk, System& sys) {
    element_init_G1(dk.D1, sys.pairing);
    for (int i = 0; i < ATTRS; i++) {
        element_init_G1(dk.Dm[i], sys.pairing);
        element_init_G1(dk.Dm_p[i], sys.pairing);
    }
}

void clear_user_dec_key(UserDecKey& dk) {
    element_clear(dk.D1);
    for (int i = 0; i < ATTRS; i++) {
        element_clear(dk.Dm[i]);
        element_clear(dk.Dm_p[i]);
    }
}

void init_cipher_block(CipherBlock& cb, System& sys) {
    element_init_Zr(cb.T1, sys.pairing);
    element_init_GT(cb.Ci, sys.pairing);
    element_init_GT(cb.c, sys.pairing);
    element_init_G1(cb.C1, sys.pairing);

    for (int i = 0; i < ATTRS; i++) {
        element_init_G1(cb.Cy[i], sys.pairing);
        element_init_G1(cb.Cy_p[i], sys.pairing);
    }
}

void clear_cipher_block(CipherBlock& cb) {
    element_clear(cb.T1);
    element_clear(cb.Ci);
    element_clear(cb.c);
    element_clear(cb.C1);

    for (int i = 0; i < ATTRS; i++) {
        element_clear(cb.Cy[i]);
        element_clear(cb.Cy_p[i]);
    }
}

void init_tag_block(TagBlock& tb, System& sys) {
    element_init_G1(tb.sigma, sys.pairing);
    for (int j = 0; j < S_TAG; j++) {
        element_init_Zr(tb.cij[j], sys.pairing);
    }
}

void clear_tag_block(TagBlock& tb) {
    element_clear(tb.sigma);
    for (int j = 0; j < S_TAG; j++) {
        element_clear(tb.cij[j]);
    }
}

void init_proof(Proof& proof, System& sys) {
    element_init_G1(proof.sigma_agg, sys.pairing);
    for (int j = 0; j < S_TAG; j++) {
        element_init_Zr(proof.mu[j], sys.pairing);
    }
}

void clear_proof(Proof& proof) {
    element_clear(proof.sigma_agg);
    for (int j = 0; j < S_TAG; j++) {
        element_clear(proof.mu[j]);
    }
}

// KeyGen: DO 和 DA 的身份签名密钥
void sign_keygen(System& sys, OwnerKey& ok, AuditorKey& ak) {
    element_t hid;

    element_init_G1(hid, sys.pairing);

    H_G1(hid, "idDO", 0, 0);
    element_pow_zn(ok.Kdos, hid, sys.s);
    element_random(ok.alphaDO);
    element_pow_zn(ok.Kdp, sys.g, ok.alphaDO);

    H_G1(hid, "idDA", 0, 0);
    element_pow_zn(ak.Kdas, hid, sys.s);
    element_random(ak.alphaDA);
    element_pow_zn(ak.Ktp, sys.g, ak.alphaDA);

    element_clear(hid);
}

// KeyGen: DU 的 CP-ABE 解密密钥
void dec_keygen(System& sys, UserDecKey& dk) {
    element_t hid, hattr, rm, inv_beta, exp, tmp;

    element_init_G1(hid, sys.pairing);
    element_init_G1(hattr, sys.pairing);
    element_init_Zr(rm, sys.pairing);
    element_init_Zr(inv_beta, sys.pairing);
    element_init_Zr(exp, sys.pairing);
    element_init_G1(tmp, sys.pairing);

    H_G1(hid, "idDU", 0, 0);

    // D1 = (g^alpha * H1(idDU))^(1/beta)
    element_invert(inv_beta, sys.beta);
    element_pow_zn(tmp, sys.g, sys.alpha);
    element_mul(tmp, tmp, hid);
    element_pow_zn(dk.D1, tmp, inv_beta);

    // Dm = H1(idDU) * H1(m)^rm, Dm' = g^rm
    for (int i = 0; i < ATTRS; i++) {
        element_random(rm);

        H_G1(hattr, "attr", i, 0);
        element_pow_zn(tmp, hattr, rm);
        element_mul(dk.Dm[i], hid, tmp);

        element_pow_zn(dk.Dm_p[i], sys.g, rm);
    }

    element_clear(hid);
    element_clear(hattr);
    element_clear(rm);
    element_clear(inv_beta);
    element_clear(exp);
    element_clear(tmp);
}

// CiphertextGen: FC = {c, C1, Cy, Cy'}
void encrypt_block(System& sys, CipherBlock& cb, int block_id) {
    element_t qx, tmpGT, tmpG, h;

    element_init_Zr(qx, sys.pairing);
    element_init_GT(tmpGT, sys.pairing);
    element_init_G1(tmpG, sys.pairing);
    element_init_G1(h, sys.pairing);

    element_random(cb.T1);
    element_random(cb.Ci);

    // c = Ci * e(g,g)^(alpha*T1)
    element_pow_zn(tmpGT, sys.egg_alpha, cb.T1);
    element_mul(cb.c, cb.Ci, tmpGT);

    // C1 = g^(beta*T1)
    element_pow_zn(cb.C1, sys.g_beta, cb.T1);

    // 叶子节点：Cy = g^q_y(0), Cy' = H1(att(y))^q_y(0)
    for (int i = 0; i < ATTRS; i++) {
        element_random(qx);
        element_pow_zn(cb.Cy[i], sys.g, qx);

        H_G1(h, "attr", i, 0);
        element_pow_zn(cb.Cy_p[i], h, qx);
    }

    element_clear(qx);
    element_clear(tmpGT);
    element_clear(tmpG);
    element_clear(h);
}

// TagGen: sigma_i = Kdos * (H4(eta||SKF||i) * product Aj^cij)^alphaDO
void tag_gen(System& sys, OwnerKey& ok, TagBlock& tb, int block_id) {
    element_t h4, base, tmp, aj;

    element_init_G1(h4, sys.pairing);
    element_init_G1(base, sys.pairing);
    element_init_G1(tmp, sys.pairing);
    element_init_G1(aj, sys.pairing);

    H_G1(h4, "H4_eta_SKF_i", block_id, 0);
    element_set(base, h4);

    for (int j = 0; j < S_TAG; j++) {
        H_Zr(tb.cij[j], "cij", block_id, j);

        // A_j = g^{a_j}，这里用 H_G1 固定模拟公共 A_j
        H_G1(aj, "Aj", j, 0);
        element_pow_zn(tmp, aj, tb.cij[j]);
        element_mul(base, base, tmp);
    }

    element_pow_zn(base, base, ok.alphaDO);
    element_mul(tb.sigma, ok.Kdos, base);

    element_clear(h4);
    element_clear(base);
    element_clear(tmp);
    element_clear(aj);
}

// ProofGen: sigma = product sigma_i^v_i, mu_j = sum v_i*cij
void proof_gen(System& sys, TagBlock& tb, Proof& proof, int block_id) {
    element_t vi, tmp;

    element_init_Zr(vi, sys.pairing);
    element_init_G1(tmp, sys.pairing);

    H_Zr(vi, "audit_v", block_id, 0);

    element_pow_zn(proof.sigma_agg, tb.sigma, vi);

    for (int j = 0; j < S_TAG; j++) {
        element_mul(proof.mu[j], vi, tb.cij[j]);
    }

    element_clear(vi);
    element_clear(tmp);
}

// ProofVerify:
// e(sigma,g) ?= e(H1(idDO)^v, MPK) * e(H4^v * product Aj^mu_j, Kdp)
void proof_verify(System& sys, OwnerKey& ok, Proof& proof, int block_id) {
    element_t left, right1, right2, right, hid, h4, vi, base, tmp, aj;

    element_init_GT(left, sys.pairing);
    element_init_GT(right1, sys.pairing);
    element_init_GT(right2, sys.pairing);
    element_init_GT(right, sys.pairing);

    element_init_G1(hid, sys.pairing);
    element_init_G1(h4, sys.pairing);
    element_init_G1(base, sys.pairing);
    element_init_G1(tmp, sys.pairing);
    element_init_G1(aj, sys.pairing);
    element_init_Zr(vi, sys.pairing);

    H_Zr(vi, "audit_v", block_id, 0);

    element_pairing(left, proof.sigma_agg, sys.g);

    H_G1(hid, "idDO", 0, 0);
    element_pow_zn(hid, hid, vi);
    element_pairing(right1, hid, sys.MPK);

    H_G1(h4, "H4_eta_SKF_i", block_id, 0);
    element_pow_zn(base, h4, vi);

    for (int j = 0; j < S_TAG; j++) {
        H_G1(aj, "Aj", j, 0);
        element_pow_zn(tmp, aj, proof.mu[j]);
        element_mul(base, base, tmp);
    }

    element_pairing(right2, base, ok.Kdp);
    element_mul(right, right1, right2);

    g_sink ^= element_cmp(left, right);

    element_clear(left);
    element_clear(right1);
    element_clear(right2);
    element_clear(right);
    element_clear(hid);
    element_clear(h4);
    element_clear(vi);
    element_clear(base);
    element_clear(tmp);
    element_clear(aj);
}

// Dec:
// Fx = e(Dm, Cy) / e(Cy', Dm')
// B = e(C1, D1)
// Ci = c * FR / B
void decrypt_block(System& sys, UserDecKey& dk, CipherBlock& cb) {
    element_t p1, p2, inv, Fx, FR, B, result;

    element_init_GT(p1, sys.pairing);
    element_init_GT(p2, sys.pairing);
    element_init_GT(inv, sys.pairing);
    element_init_GT(Fx, sys.pairing);
    element_init_GT(FR, sys.pairing);
    element_init_GT(B, sys.pairing);
    element_init_GT(result, sys.pairing);

    element_set1(FR);

    for (int i = 0; i < ATTRS; i++) {
        element_pairing(p1, dk.Dm[i], cb.Cy[i]);
        element_pairing(p2, cb.Cy_p[i], dk.Dm_p[i]);
        element_invert(inv, p2);
        element_mul(Fx, p1, inv);
        element_mul(FR, FR, Fx);
    }

    element_pairing(B, cb.C1, dk.D1);
    element_invert(inv, B);
    element_mul(result, cb.c, FR);
    element_mul(result, result, inv);

    g_sink ^= element_is1(result);

    element_clear(p1);
    element_clear(p2);
    element_clear(inv);
    element_clear(Fx);
    element_clear(FR);
    element_clear(B);
    element_clear(result);
}

int main() {
    std::printf("blocks,dec_keygen_ms,sign_keygen_ms,encryption_ms,signature_ms,decryption_ms,verification_ms\n");

    for (int blocks = 50; blocks <= 1000; blocks += 50) {
        System sys;
        init_system(sys);

        OwnerKey ok;
        AuditorKey ak;
        UserDecKey dk;

        init_owner_key(ok, sys);
        init_auditor_key(ak, sys);
        init_user_dec_key(dk, sys);

        auto st = Clock::now();
        dec_keygen(sys, dk);
        auto ed = Clock::now();
        double dec_keygen_ms = ms_since(st, ed);

        st = Clock::now();
        sign_keygen(sys, ok, ak);
        ed = Clock::now();
        double sign_keygen_ms = ms_since(st, ed);

        double encryption_ms = 0.0;
        double signature_ms = 0.0;
        double decryption_ms = 0.0;
        double verification_ms = 0.0;

        for (int b = 0; b < blocks; b++) {
            CipherBlock cb;
            TagBlock tb;
            Proof proof;

            init_cipher_block(cb, sys);
            init_tag_block(tb, sys);
            init_proof(proof, sys);

            st = Clock::now();
            encrypt_block(sys, cb, b);
            ed = Clock::now();
            encryption_ms += ms_since(st, ed);

            st = Clock::now();
            tag_gen(sys, ok, tb, b);
            ed = Clock::now();
            signature_ms += ms_since(st, ed);

            st = Clock::now();
            decrypt_block(sys, dk, cb);
            ed = Clock::now();
            decryption_ms += ms_since(st, ed);

            st = Clock::now();
            proof_gen(sys, tb, proof, b);
            proof_verify(sys, ok, proof, b);
            ed = Clock::now();
            verification_ms += ms_since(st, ed);

            clear_proof(proof);
            clear_tag_block(tb);
            clear_cipher_block(cb);
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

        clear_user_dec_key(dk);
        clear_auditor_key(ak);
        clear_owner_key(ok);
        clear_system(sys);
    }

    if (g_sink == 123456789) {
        std::printf("ignore: %d\n", g_sink);
    }

    return 0;
}