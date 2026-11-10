#define _CRT_SECURE_NO_WARNINGS
#include <pbc.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using Clock = std::chrono::high_resolution_clock;

// 固定参与用户数量为 1
static const int PARTICIPANTS = 1;

// Wang 论文是多方签名方案，不是加密方案
// 为了和前面实验保持同样输出格式，仍然输出六列时间

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

// 模拟哈希到 Zr
void H_Zr(element_t out, const char* tag, int block_id, int extra) {
    char buf[128];
    std::snprintf(buf, sizeof(buf), "%s_%d_%d", tag, block_id, extra);
    element_from_hash(out, buf, (int)std::strlen(buf));
}

struct System {
    pairing_t pairing;

    // G 群生成元
    element_t g;

    // 承诺参数，模拟论文中的乘法同态承诺 key
    element_t h;
};

struct UserKey {
    // 私钥片段 x_i
    element_t sk;

    // 公钥片段 X_i = g^{x_i}
    element_t pk;
};

struct Commit {
    // 第一轮承诺 y_i = g^{r_i}
    element_t y;

    // 承诺随机数 c_i，这里用 G1 元素模拟
    element_t c;

    // 签名随机数 r_i
    element_t r;
};

struct CRFCommit {
    // CRF 输入侧重随机化后的承诺 y_i'
    element_t y_prime;

    // CRF 输入侧重随机化后的承诺 c_i'
    element_t c_prime;

    // CRF 使用的随机量 rho_i
    element_t rho;
};

struct Signature {
    // 聚合承诺 R
    element_t R;

    // 聚合公钥 PK
    element_t APK;

    // challenge e = H(m, R, APK)
    element_t e;

    // 签名响应 z
    element_t z;

    // CRF 输出侧重随机化量
    element_t eta;
};

void init_system(System& sys) {
    pairing_init_set_buf(sys.pairing, DEFAULT_PARAM, std::strlen(DEFAULT_PARAM));

    element_init_G1(sys.g, sys.pairing);
    element_init_G1(sys.h, sys.pairing);

    element_random(sys.g);
    element_random(sys.h);
}

void clear_system(System& sys) {
    element_clear(sys.g);
    element_clear(sys.h);
    pairing_clear(sys.pairing);
}

void init_user_key(UserKey& key, System& sys) {
    element_init_Zr(key.sk, sys.pairing);
    element_init_G1(key.pk, sys.pairing);
}

void clear_user_key(UserKey& key) {
    element_clear(key.sk);
    element_clear(key.pk);
}

void init_commit(Commit& cm, System& sys) {
    element_init_G1(cm.y, sys.pairing);
    element_init_G1(cm.c, sys.pairing);
    element_init_Zr(cm.r, sys.pairing);
}

void clear_commit(Commit& cm) {
    element_clear(cm.y);
    element_clear(cm.c);
    element_clear(cm.r);
}

void init_crf_commit(CRFCommit& cc, System& sys) {
    element_init_G1(cc.y_prime, sys.pairing);
    element_init_G1(cc.c_prime, sys.pairing);
    element_init_Zr(cc.rho, sys.pairing);
}

void clear_crf_commit(CRFCommit& cc) {
    element_clear(cc.y_prime);
    element_clear(cc.c_prime);
    element_clear(cc.rho);
}

void init_signature(Signature& sig, System& sys) {
    element_init_G1(sig.R, sys.pairing);
    element_init_G1(sig.APK, sys.pairing);
    element_init_Zr(sig.e, sys.pairing);
    element_init_Zr(sig.z, sys.pairing);
    element_init_Zr(sig.eta, sys.pairing);
}

void clear_signature(Signature& sig) {
    element_clear(sig.R);
    element_clear(sig.APK);
    element_clear(sig.e);
    element_clear(sig.z);
    element_clear(sig.eta);
}

// KeyGen：每个参与者生成 sk_i 和 pk_i = g^{sk_i}
void keygen(System& sys, UserKey* keys) {
    for (int i = 0; i < PARTICIPANTS; i++) {
        element_random(keys[i].sk);
        element_pow_zn(keys[i].pk, sys.g, keys[i].sk);
    }
}

// Phase 1：参与者生成承诺 y_i = g^{r_i}, c_i = h^{r_i}
void phase1_commit(System& sys, Commit* cm, int block_id) {
    for (int i = 0; i < PARTICIPANTS; i++) {
        element_random(cm[i].r);

        element_pow_zn(cm[i].y, sys.g, cm[i].r);
        element_pow_zn(cm[i].c, sys.h, cm[i].r);

        // 加入一次 hash，模拟消息绑定
        element_t msg_hash;
        element_init_Zr(msg_hash, sys.pairing);
        H_Zr(msg_hash, "MCRF_MESSAGE", block_id, i);
        element_pow_zn(cm[i].c, cm[i].c, msg_hash);
        element_clear(msg_hash);
    }
}

// CRF 输入侧处理：y_i' = y_i * g^{rho_i}, c_i' = c_i * h^{rho_i}
void input_crf(System& sys, Commit* cm, CRFCommit* cc) {
    element_t tmp;

    element_init_G1(tmp, sys.pairing);

    for (int i = 0; i < PARTICIPANTS; i++) {
        element_random(cc[i].rho);

        element_pow_zn(tmp, sys.g, cc[i].rho);
        element_mul(cc[i].y_prime, cm[i].y, tmp);

        element_pow_zn(tmp, sys.h, cc[i].rho);
        element_mul(cc[i].c_prime, cm[i].c, tmp);
    }

    element_clear(tmp);
}

// Phase 2/3：聚合承诺和公钥，计算 challenge e
void aggregate_and_challenge(System& sys, UserKey* keys, CRFCommit* cc, Signature& sig, int block_id) {
    element_set1(sig.R);
    element_set1(sig.APK);

    for (int i = 0; i < PARTICIPANTS; i++) {
        element_mul(sig.R, sig.R, cc[i].y_prime);
        element_mul(sig.APK, sig.APK, keys[i].pk);
    }

    H_Zr(sig.e, "MCRF_CHALLENGE", block_id, 0);
}

// Phase 4：计算签名响应 z_i = r_i + rho_i + e * sk_i
// 聚合得到 z = sum z_i
void response_phase(System& sys, UserKey* keys, Commit* cm, CRFCommit* cc, Signature& sig) {
    element_t zi, tmp;

    element_init_Zr(zi, sys.pairing);
    element_init_Zr(tmp, sys.pairing);

    element_set0(sig.z);

    for (int i = 0; i < PARTICIPANTS; i++) {
        element_mul(tmp, sig.e, keys[i].sk);
        element_add(zi, cm[i].r, cc[i].rho);
        element_add(zi, zi, tmp);

        element_add(sig.z, sig.z, zi);
    }

    element_clear(zi);
    element_clear(tmp);
}

// CRF 输出侧处理：对最终响应再次重随机化
// 这里模拟论文中的 output-side CRF re-randomization
void output_crf(System& sys, Signature& sig) {
    element_random(sig.eta);
    element_add(sig.z, sig.z, sig.eta);

    element_t tmp;
    element_init_G1(tmp, sys.pairing);

    // R' = R * g^{eta}，保持验证等式成立
    element_pow_zn(tmp, sys.g, sig.eta);
    element_mul(sig.R, sig.R, tmp);

    element_clear(tmp);
}

// MSign：完整多方签名过程
void mcrf_sign(System& sys, UserKey* keys, Commit* cm, CRFCommit* cc, Signature& sig, int block_id) {
    phase1_commit(sys, cm, block_id);
    input_crf(sys, cm, cc);
    aggregate_and_challenge(sys, keys, cc, sig, block_id);
    response_phase(sys, keys, cm, cc, sig);
    output_crf(sys, sig);
}

// Verify：检查 g^z ?= R * APK^e
int verify(System& sys, Signature& sig) {
    element_t left, right, tmp;

    element_init_G1(left, sys.pairing);
    element_init_G1(right, sys.pairing);
    element_init_G1(tmp, sys.pairing);

    element_pow_zn(left, sys.g, sig.z);
    element_pow_zn(tmp, sig.APK, sig.e);
    element_mul(right, sig.R, tmp);

    int ok = (element_cmp(left, right) == 0);

    g_sink ^= ok;

    element_clear(left);
    element_clear(right);
    element_clear(tmp);

    return ok;
}

int main() {
    std::printf("blocks,dec_keygen_ms,sign_keygen_ms,encryption_ms,signature_ms,decryption_ms,verification_ms\n");

    for (int blocks = 50; blocks <= 1000; blocks += 50) {
        System sys;
        init_system(sys);

        UserKey keys[PARTICIPANTS];
        for (int i = 0; i < PARTICIPANTS; i++) {
            init_user_key(keys[i], sys);
        }

        // Wang 论文没有解密密钥生成
        double dec_keygen_ms = 0.0;

        auto st = Clock::now();
        keygen(sys, keys);
        auto ed = Clock::now();
        double sign_keygen_ms = ms_since(st, ed);

        // Wang 论文没有加密、解密算法
        double encryption_ms = 0.0;
        double decryption_ms = 0.0;

        double signature_ms = 0.0;
        double verification_ms = 0.0;

        for (int b = 0; b < blocks; b++) {
            Commit cm[PARTICIPANTS];
            CRFCommit cc[PARTICIPANTS];
            Signature sig;

            for (int i = 0; i < PARTICIPANTS; i++) {
                init_commit(cm[i], sys);
                init_crf_commit(cc[i], sys);
            }
            init_signature(sig, sys);

            st = Clock::now();
            mcrf_sign(sys, keys, cm, cc, sig, b);
            ed = Clock::now();
            signature_ms += ms_since(st, ed);

            st = Clock::now();
            verify(sys, sig);
            ed = Clock::now();
            verification_ms += ms_since(st, ed);

            clear_signature(sig);
            for (int i = 0; i < PARTICIPANTS; i++) {
                clear_crf_commit(cc[i]);
                clear_commit(cm[i]);
            }
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

        for (int i = 0; i < PARTICIPANTS; i++) {
            clear_user_key(keys[i]);
        }
        clear_system(sys);
    }

    if (g_sink == 123456789) {
        std::printf("ignore: %d\n", g_sink);
    }

    return 0;
}