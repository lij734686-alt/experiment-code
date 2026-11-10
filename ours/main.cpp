
#include <pbc.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>

/*
 * - 属性全集 U
 * - 签名者属性集合 W
 * - 审计者属性集合 W_A
 * - LSSS 矩阵 M 和 rho 映射
 * - 审计者解密密钥 ASK
 * - 文件授权密文 CT_j
 * - 文件绑定属性签名密钥 SSK
 * - 每个数据块签名 sigma_{j,t}
 * - 授权解密恢复值 Z_j
 * - 审计响应和验证侧中间元素
 *

 */

struct SchemeParams {
    int n_attrs;      // n：属性数量
    int l_rows;       // l：LSSS 矩阵行数
    int poly_degree;  // L：R_{j,t}=P0 prod P_l^{H3(...)} 中的 L
    int file_count;   // N：满足关键词和策略的文件数量
};

struct ZrVector {
    element_t* items;
    int count;
};

struct LsssPolicy {
    element_t** m;
    int rows;
    int cols;
    int* rho;
};

struct PublicContext {
    pairing_ptr pairing;
    SchemeParams params;

    element_t g;
    element_t u;
    element_t alpha_s;
    element_t alpha_A;
    element_t beta;
    element_t g_beta;
    element_t pk_t;
    element_t egg_alpha_A;
    element_t h0;

    element_t* h_attrs;
    element_t* p_list;
};

struct AuditorKey {
    element_t K_A;
    element_t L_A;
    element_t* K_A_x;
    int* attrs;
    int attr_count;
};

struct FileCiphertext {
    element_t fid;
    element_t x_j;
    element_t X_j;
    element_t theta_j;
    element_t Y_j;
    element_t HF_j;
    element_t h0_j;
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

struct AuditResponse {
    element_t* M_j;
    element_t* delta1_j;
    element_t delta2;
    element_t** sigma0;
};

using Clock = std::chrono::high_resolution_clock;

static double elapsed_ms(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

static void die_alloc(void* ptr) {
    if (!ptr) {
        std::fprintf(stderr, "内存分配失败。\n");
        std::exit(1);
    }
}

static element_t* alloc_elements(int count) {
    element_t* items = static_cast<element_t*>(std::calloc(static_cast<size_t>(count), sizeof(element_t)));
    die_alloc(items);
    return items;
}

static void init_zr_vector(ZrVector* vec, pairing_ptr pairing, int count, element_t first) {
    vec->count = count;
    vec->items = alloc_elements(count);
    for (int i = 0; i < count; i++) {
        element_init_Zr(vec->items[i], pairing);
        if (i == 0) {
            element_set(vec->items[i], first);
        }
        else {
            element_random(vec->items[i]);
        }
    }
}

static void clear_zr_vector(ZrVector* vec) {
    for (int i = 0; i < vec->count; i++) {
        element_clear(vec->items[i]);
    }
    std::free(vec->items);
}

static int* make_full_attribute_set(int n_attrs) {
    int* attrs = static_cast<int*>(std::calloc(static_cast<size_t>(n_attrs), sizeof(int)));
    die_alloc(attrs);
    for (int i = 0; i < n_attrs; i++) {
        attrs[i] = i;
    }
    return attrs;
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
    element_init_GT(ctx->pk_t, pairing);
    element_init_GT(ctx->egg_alpha_A, pairing);
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
    element_pow_zn(ctx->pk_t, egg, ctx->alpha_s);
    element_pow_zn(ctx->egg_alpha_A, egg, ctx->alpha_A);
    element_clear(egg);

    ctx->h_attrs = alloc_elements(params.n_attrs);
    for (int i = 0; i < params.n_attrs; i++) {
        element_init_G1(ctx->h_attrs[i], pairing);
        element_random(ctx->h_attrs[i]);
    }

    ctx->p_list = alloc_elements(params.poly_degree + 1);
    for (int i = 0; i <= params.poly_degree; i++) {
        element_t x_l;
        element_init_Zr(x_l, pairing);
        element_random(x_l);
        element_init_G1(ctx->p_list[i], pairing);
        element_pow_zn(ctx->p_list[i], ctx->g, x_l);
        element_clear(x_l);
    }
}

static void clear_public_context(PublicContext* ctx) {
    element_clear(ctx->g);
    element_clear(ctx->u);
    element_clear(ctx->alpha_s);
    element_clear(ctx->alpha_A);
    element_clear(ctx->beta);
    element_clear(ctx->g_beta);
    element_clear(ctx->pk_t);
    element_clear(ctx->egg_alpha_A);
    element_clear(ctx->h0);

    for (int i = 0; i < ctx->params.n_attrs; i++) {
        element_clear(ctx->h_attrs[i]);
    }
    std::free(ctx->h_attrs);

    for (int i = 0; i <= ctx->params.poly_degree; i++) {
        element_clear(ctx->p_list[i]);
    }
    std::free(ctx->p_list);
}

static void init_policy(LsssPolicy* policy, pairing_ptr pairing, int rows, int cols, int n_attrs) {
    policy->rows = rows;
    policy->cols = cols;
    policy->rho = static_cast<int*>(std::calloc(static_cast<size_t>(rows), sizeof(int)));
    policy->m = static_cast<element_t**>(std::calloc(static_cast<size_t>(rows), sizeof(element_t*)));
    die_alloc(policy->rho);
    die_alloc(policy->m);

    for (int i = 0; i < rows; i++) {
        policy->m[i] = alloc_elements(cols);
        policy->rho[i] = i % n_attrs;
        for (int j = 0; j < cols; j++) {
            element_init_Zr(policy->m[i][j], pairing);
            element_set0(policy->m[i][j]);
        }
    }

    // 构造一个可满足的 LSSS：第一行是 (1,0,...,0)，其他行用单位向量循环补齐。
    element_set1(policy->m[0][0]);
    for (int i = 1; i < rows; i++) {
        element_set1(policy->m[i][i % cols]);
    }
}

static void clear_policy(LsssPolicy* policy) {
    for (int i = 0; i < policy->rows; i++) {
        for (int j = 0; j < policy->cols; j++) {
            element_clear(policy->m[i][j]);
        }
        std::free(policy->m[i]);
    }
    std::free(policy->m);
    std::free(policy->rho);
}

static void dot_row_vector(element_t out, pairing_ptr pairing, LsssPolicy* policy, int row, ZrVector* vec) {
    element_set0(out);
    for (int j = 0; j < policy->cols; j++) {
        element_t tmp;
        element_init_Zr(tmp, pairing);
        element_mul(tmp, policy->m[row][j], vec->items[j]);
        element_add(out, out, tmp);
        element_clear(tmp);
    }
}

static void init_auditor_key(AuditorKey* key, PublicContext* ctx) {
    int n = ctx->params.n_attrs;
    pairing_ptr pairing = ctx->pairing;

    key->attr_count = n;
    key->attrs = make_full_attribute_set(n);
    key->K_A_x = alloc_elements(n);

    element_init_G1(key->K_A, pairing);
    element_init_G1(key->L_A, pairing);

    element_t t_A, part1, part2;
    element_init_Zr(t_A, pairing);
    element_init_G1(part1, pairing);
    element_init_G1(part2, pairing);
    element_random(t_A);

    // K_A = g^{alpha_A}(g^beta)^{t_A}
    element_pow_zn(part1, ctx->g, ctx->alpha_A);
    element_pow_zn(part2, ctx->g_beta, t_A);
    element_mul(key->K_A, part1, part2);

    // L_A = g^{t_A}
    element_pow_zn(key->L_A, ctx->g, t_A);

    // K_{A,x}=h_x^{t_A}
    for (int i = 0; i < n; i++) {
        element_init_G1(key->K_A_x[i], pairing);
        element_pow_zn(key->K_A_x[i], ctx->h_attrs[i], t_A);
    }

    element_clear(t_A);
    element_clear(part1);
    element_clear(part2);
}

static void clear_auditor_key(AuditorKey* key) {
    element_clear(key->K_A);
    element_clear(key->L_A);
    for (int i = 0; i < key->attr_count; i++) {
        element_clear(key->K_A_x[i]);
    }
    std::free(key->K_A_x);
    std::free(key->attrs);
}

static void init_file_ciphertext(FileCiphertext* ct, PublicContext* ctx, LsssPolicy* policy) {
    int l = policy->rows;
    pairing_ptr pairing = ctx->pairing;

    element_init_Zr(ct->fid, pairing);
    element_init_Zr(ct->x_j, pairing);
    element_init_G1(ct->X_j, pairing);
    element_init_Zr(ct->theta_j, pairing);
    element_init_G1(ct->Y_j, pairing);
    element_init_G1(ct->HF_j, pairing);
    element_init_G1(ct->h0_j, pairing);
    element_init_Zr(ct->s_j, pairing);
    element_init_G1(ct->C_hat_prime, pairing);

    element_random(ct->fid);
    element_random(ct->x_j);
    element_pow_zn(ct->X_j, ctx->g, ct->x_j);
    element_random(ct->theta_j);
    element_pow_zn(ct->Y_j, ctx->g, ct->theta_j);
    element_random(ct->HF_j);
    element_random(ct->h0_j);
    element_random(ct->s_j);
    element_pow_zn(ct->C_hat_prime, ctx->g, ct->s_j);

    ZrVector y_vec;
    init_zr_vector(&y_vec, pairing, policy->cols, ct->s_j);

    ct->C_i = alloc_elements(l);
    ct->D_i = alloc_elements(l);

    for (int i = 0; i < l; i++) {
        element_t lambda_i, r_i, left, right, right_inv;
        element_init_Zr(lambda_i, pairing);
        element_init_Zr(r_i, pairing);
        element_init_G1(left, pairing);
        element_init_G1(right, pairing);
        element_init_G1(right_inv, pairing);
        element_random(r_i);

        dot_row_vector(lambda_i, pairing, policy, i, &y_vec);

        // C_i=(g^beta)^{lambda_i} h_{rho(i)}^{-r_i}
        element_init_G1(ct->C_i[i], pairing);
        element_init_G1(ct->D_i[i], pairing);
        element_pow_zn(left, ctx->g_beta, lambda_i);
        element_pow_zn(right, ctx->h_attrs[policy->rho[i]], r_i);
        element_invert(right_inv, right);
        element_mul(ct->C_i[i], left, right_inv);

        // D_i=g^{r_i}
        element_pow_zn(ct->D_i[i], ctx->g, r_i);

        element_clear(lambda_i);
        element_clear(r_i);
        element_clear(left);
        element_clear(right);
        element_clear(right_inv);
    }

    clear_zr_vector(&y_vec);
}

static void clear_file_ciphertext(FileCiphertext* ct, LsssPolicy* policy) {
    element_clear(ct->fid);
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

static void init_signing_key(SigningKey* sk, PublicContext* ctx, LsssPolicy* policy, FileCiphertext* ct) {
    int l = policy->rows;
    int n = ctx->params.n_attrs;
    pairing_ptr pairing = ctx->pairing;

    ZrVector v_vec;
    init_zr_vector(&v_vec, pairing, policy->cols, ctx->alpha_s);

    sk->K_i = alloc_elements(l);
    sk->K_i_prime = alloc_elements(l);
    sk->K_i_x_double_prime = static_cast<element_t**>(std::calloc(static_cast<size_t>(l), sizeof(element_t*)));
    die_alloc(sk->K_i_x_double_prime);

    for (int i = 0; i < l; i++) {
        element_t lambda_i, u_i, base, part1, part2;
        element_init_Zr(lambda_i, pairing);
        element_init_Zr(u_i, pairing);
        element_init_G1(base, pairing);
        element_init_G1(part1, pairing);
        element_init_G1(part2, pairing);
        element_random(u_i);

        dot_row_vector(lambda_i, pairing, policy, i, &v_vec);

        // K_i=g^{lambda_i}(h_{0,j}h_{rho(i)})^{u_i}
        element_init_G1(sk->K_i[i], pairing);
        element_mul(base, ct->h0_j, ctx->h_attrs[policy->rho[i]]);
        element_pow_zn(part1, ctx->g, lambda_i);
        element_pow_zn(part2, base, u_i);
        element_mul(sk->K_i[i], part1, part2);

        // K_i'=g^{u_i}
        element_init_G1(sk->K_i_prime[i], pairing);
        element_pow_zn(sk->K_i_prime[i], ctx->g, u_i);

        // K''_{i,x}=h_x^{u_i}, x != rho(i)
        sk->K_i_x_double_prime[i] = alloc_elements(n);
        for (int x = 0; x < n; x++) {
            element_init_G1(sk->K_i_x_double_prime[i][x], pairing);
            if (x == policy->rho[i]) {
                element_set1(sk->K_i_x_double_prime[i][x]);
            }
            else {
                element_pow_zn(sk->K_i_x_double_prime[i][x], ctx->h_attrs[x], u_i);
            }
        }

        element_clear(lambda_i);
        element_clear(u_i);
        element_clear(base);
        element_clear(part1);
        element_clear(part2);
    }

    clear_zr_vector(&v_vec);
}

static void clear_signing_key(SigningKey* sk, LsssPolicy* policy, int n_attrs) {
    for (int i = 0; i < policy->rows; i++) {
        element_clear(sk->K_i[i]);
        element_clear(sk->K_i_prime[i]);
        for (int x = 0; x < n_attrs; x++) {
            element_clear(sk->K_i_x_double_prime[i][x]);
        }
        std::free(sk->K_i_x_double_prime[i]);
    }
    std::free(sk->K_i);
    std::free(sk->K_i_prime);
    std::free(sk->K_i_x_double_prime);
}

static void compute_R(element_t out, PublicContext* ctx) {
    element_set(out, ctx->p_list[0]);
    for (int ell = 1; ell <= ctx->params.poly_degree; ell++) {
        element_t h, part;
        element_init_Zr(h, ctx->pairing);
        element_init_G1(part, ctx->pairing);
        element_random(h);
        element_pow_zn(part, ctx->p_list[ell], h);
        element_mul(out, out, part);
        element_clear(h);
        element_clear(part);
    }
}

static void init_block_signatures(BlockSignature* sigs, int block_count,
    PublicContext* ctx, LsssPolicy* policy,
    FileCiphertext* ct, SigningKey* sk) {
    int n = ctx->params.n_attrs;
    pairing_ptr pairing = ctx->pairing;

    for (int t = 0; t < block_count; t++) {
        element_t xi, r, R, HW, tmp, product, msg_part, r_part;
        element_init_Zr(xi, pairing);
        element_init_Zr(r, pairing);
        element_init_G1(R, pairing);
        element_init_G1(HW, pairing);
        element_init_G1(tmp, pairing);
        element_init_G1(product, pairing);
        element_init_G1(msg_part, pairing);
        element_init_G1(r_part, pairing);

        element_random(xi);
        element_random(r);
        compute_R(R, ctx);

        element_init_G1(sigs[t].sigma0, pairing);
        element_init_G1(sigs[t].sigma1, pairing);
        element_init_G1(sigs[t].sigma2, pairing);

        // sigma0=g^r
        element_pow_zn(sigs[t].sigma0, ctx->g, r);

        // 示例策略中 I={0}, omega_0=1，所以 sigma1=K'_0 * g^xi。
        element_pow_zn(tmp, ctx->g, xi);
        element_mul(sigs[t].sigma1, sk->K_i_prime[0], tmp);

        // HW_j=h_{0,j} prod_{att_x in W} h_x
        element_set(HW, ct->h0_j);
        for (int x = 0; x < n; x++) {
            element_mul(HW, HW, ctx->h_attrs[x]);
        }

        // sigma2 的主要群元素构造。
        element_set(product, sk->K_i[0]);
        for (int x = 0; x < n; x++) {
            if (x != policy->rho[0]) {
                element_mul(product, product, sk->K_i_x_double_prime[0][x]);
            }
        }

        element_pow_zn(tmp, HW, xi);
        element_mul(product, product, tmp);
        element_mul(product, product, ct->HF_j);

        // 模拟 (H1(ctx) * u^c)^{x_j}
        element_random(msg_part);
        element_pow_zn(msg_part, msg_part, ct->x_j);
        element_mul(product, product, msg_part);

        // R_{j,t}^{r_{j,t}}
        element_pow_zn(r_part, R, r);
        element_mul(sigs[t].sigma2, product, r_part);

        element_clear(xi);
        element_clear(r);
        element_clear(R);
        element_clear(HW);
        element_clear(tmp);
        element_clear(product);
        element_clear(msg_part);
        element_clear(r_part);
    }
}

static void clear_block_signatures(BlockSignature* sigs, int block_count) {
    for (int t = 0; t < block_count; t++) {
        element_clear(sigs[t].sigma0);
        element_clear(sigs[t].sigma1);
        element_clear(sigs[t].sigma2);
    }
}

static void authorization_decrypt(element_t out, PublicContext* ctx,
    LsssPolicy* policy, AuditorKey* ask,
    FileCiphertext* ct) {
    pairing_ptr pairing = ctx->pairing;
    element_t numerator, denominator, pair1, pair2, tmp;

    element_init_GT(numerator, pairing);
    element_init_GT(denominator, pairing);
    element_init_GT(pair1, pairing);
    element_init_GT(pair2, pairing);
    element_init_GT(tmp, pairing);
    element_set1(denominator);

    pairing_apply(numerator, ct->C_hat_prime, ask->K_A, pairing);

    // 属性集合包含全部属性，所以这里遍历所有 LSSS 行。
    for (int i = 0; i < policy->rows; i++) {
        pairing_apply(pair1, ct->C_i[i], ask->L_A, pairing);
        pairing_apply(pair2, ct->D_i[i], ask->K_A_x[policy->rho[i]], pairing);
        element_mul(tmp, pair1, pair2);
        element_mul(denominator, denominator, tmp);
    }

    element_invert(denominator, denominator);
    element_mul(out, numerator, denominator);

    element_clear(numerator);
    element_clear(denominator);
    element_clear(pair1);
    element_clear(pair2);
    element_clear(tmp);
}

static void init_audit_response(AuditResponse* resp, int block_count,
    PublicContext* ctx, BlockSignature** sigs) {
    int file_count = ctx->params.file_count;
    pairing_ptr pairing = ctx->pairing;

    resp->M_j = alloc_elements(file_count);
    resp->delta1_j = alloc_elements(file_count);
    element_init_G1(resp->delta2, pairing);
    element_set1(resp->delta2);

    resp->sigma0 = static_cast<element_t**>(std::calloc(static_cast<size_t>(file_count), sizeof(element_t*)));
    die_alloc(resp->sigma0);

    for (int j = 0; j < file_count; j++) {
        element_init_Zr(resp->M_j[j], pairing);
        element_random(resp->M_j[j]);
        element_init_G1(resp->delta1_j[j], pairing);
        element_set1(resp->delta1_j[j]);

        resp->sigma0[j] = alloc_elements(block_count);
        for (int t = 0; t < block_count; t++) {
            element_t v_t, part1, part2;
            element_init_Zr(v_t, pairing);
            element_init_G1(part1, pairing);
            element_init_G1(part2, pairing);
            element_random(v_t);

            element_pow_zn(part1, sigs[j][t].sigma1, v_t);
            element_mul(resp->delta1_j[j], resp->delta1_j[j], part1);

            element_pow_zn(part2, sigs[j][t].sigma2, v_t);
            element_mul(resp->delta2, resp->delta2, part2);

            element_init_G1(resp->sigma0[j][t], pairing);
            element_set(resp->sigma0[j][t], sigs[j][t].sigma0);

            element_clear(v_t);
            element_clear(part1);
            element_clear(part2);
        }
    }
}

static void clear_audit_response(AuditResponse* resp, int file_count, int block_count) {
    for (int j = 0; j < file_count; j++) {
        element_clear(resp->M_j[j]);
        element_clear(resp->delta1_j[j]);
        for (int t = 0; t < block_count; t++) {
            element_clear(resp->sigma0[j][t]);
        }
        std::free(resp->sigma0[j]);
    }
    std::free(resp->M_j);
    std::free(resp->delta1_j);
    std::free(resp->sigma0);
    element_clear(resp->delta2);
}

static int audit_verify(PublicContext* ctx, FileCiphertext* cts,
    BlockSignature** sigs, AuditResponse* resp,
    int block_count) {
    pairing_ptr pairing = ctx->pairing;
    int file_count = ctx->params.file_count;
    element_t left, right, tmp_gt, tmp_g1, exp;

    element_init_GT(left, pairing);
    element_init_GT(right, pairing);
    element_init_GT(tmp_gt, pairing);
    element_init_G1(tmp_g1, pairing);
    element_init_Zr(exp, pairing);

    pairing_apply(left, resp->delta2, ctx->g, pairing);
    element_set1(right);

    // (PK_T)^{V|S|}，这里用随机 exp 模拟 V|S|。
    element_random(exp);
    element_pow_zn(tmp_gt, ctx->pk_t, exp);
    element_mul(right, right, tmp_gt);

    for (int j = 0; j < file_count; j++) {
        element_t HW;
        element_init_G1(HW, pairing);
        element_set(HW, cts[j].h0_j);
        for (int x = 0; x < ctx->params.n_attrs; x++) {
            element_mul(HW, HW, ctx->h_attrs[x]);
        }

        pairing_apply(tmp_gt, HW, resp->delta1_j[j], pairing);
        element_mul(right, right, tmp_gt);

        element_random(exp);
        element_pow_zn(tmp_g1, cts[j].HF_j, exp);
        pairing_apply(tmp_gt, tmp_g1, ctx->g, pairing);
        element_mul(right, right, tmp_gt);

        element_random(tmp_g1);
        pairing_apply(tmp_gt, tmp_g1, cts[j].X_j, pairing);
        element_mul(right, right, tmp_gt);

        for (int t = 0; t < block_count; t++) {
            element_t R;
            element_init_G1(R, pairing);
            compute_R(R, ctx);
            pairing_apply(tmp_gt, R, sigs[j][t].sigma0, pairing);
            element_mul(right, right, tmp_gt);
            element_clear(R);
        }

        element_clear(HW);
    }

    int ok = !element_cmp(left, right);

    element_clear(left);
    element_clear(right);
    element_clear(tmp_gt);
    element_clear(tmp_g1);
    element_clear(exp);
    return ok;
}

static double time_dec_keygen(PublicContext* ctx, AuditorKey* out) {
    auto start = Clock::now();
    init_auditor_key(out, ctx);
    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static double time_encrypt(PublicContext* ctx, LsssPolicy* policy, FileCiphertext* cts) {
    auto start = Clock::now();
    for (int j = 0; j < ctx->params.file_count; j++) {
        init_file_ciphertext(&cts[j], ctx, policy);
    }
    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static double time_sign_keygen(PublicContext* ctx, LsssPolicy* policy, FileCiphertext* cts, SigningKey* sks) {
    auto start = Clock::now();
    for (int j = 0; j < ctx->params.file_count; j++) {
        init_signing_key(&sks[j], ctx, policy, &cts[j]);
    }
    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static double time_signature(PublicContext* ctx, LsssPolicy* policy,
    FileCiphertext* cts, SigningKey* sks,
    BlockSignature** sigs, int block_count) {
    auto start = Clock::now();
    for (int j = 0; j < ctx->params.file_count; j++) {
        init_block_signatures(sigs[j], block_count, ctx, policy, &cts[j], &sks[j]);
    }
    auto end = Clock::now();
    return elapsed_ms(start, end);
}

static double time_decryption(PublicContext* ctx, LsssPolicy* policy, AuditorKey* ask, FileCiphertext* cts) {
    element_t z;
    element_init_GT(z, ctx->pairing);

    auto start = Clock::now();
    for (int j = 0; j < ctx->params.file_count; j++) {
        authorization_decrypt(z, ctx, policy, ask, &cts[j]);
    }
    auto end = Clock::now();

    element_clear(z);
    return elapsed_ms(start, end);
}

static double time_verification(PublicContext* ctx, FileCiphertext* cts,
    BlockSignature** sigs, int block_count) {
    AuditResponse resp;
    init_audit_response(&resp, block_count, ctx, sigs);

    auto start = Clock::now();
    (void)audit_verify(ctx, cts, sigs, &resp, block_count);
    auto end = Clock::now();

    clear_audit_response(&resp, ctx->params.file_count, block_count);
    return elapsed_ms(start, end);
}

static void usage(const char* name) {
    std::fprintf(stderr,
        "Usage: %s <pairing.param> [n_attrs] [l_rows] [poly_degree] [file_count]\n"
        "Example: %s a.param 10 10 5 1\n",
        name, name);
}

int main(int argc, char** argv) {
    pairing_t pairing;
    char param[8192];
    FILE* fp;
    SchemeParams params = { 10, 10, 5, 1 };
    int block_sizes[] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
    int block_size_count = static_cast<int>(sizeof(block_sizes) / sizeof(block_sizes[0]));

    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    if (argc > 2) params.n_attrs = std::atoi(argv[2]);
    if (argc > 3) params.l_rows = std::atoi(argv[3]);
    if (argc > 4) params.poly_degree = std::atoi(argv[4]);
    if (argc > 5) params.file_count = std::atoi(argv[5]);

    if (params.n_attrs <= 0 || params.l_rows <= 0 || params.poly_degree < 0 || params.file_count <= 0) {
        std::fprintf(stderr, "参数必须满足 n_attrs>0, l_rows>0, poly_degree>=0, file_count>0。\n");
        return 1;
    }

    errno_t err = fopen_s(&fp, argv[1], "r");
    if (err != 0 || !fp) {
        std::perror("无法打开 PBC 参数文件");
        return 1;
    }

    size_t count = std::fread(param, 1, sizeof(param), fp);
    std::fclose(fp);
    if (!count) {
        std::fprintf(stderr, "PBC 参数文件为空。\n");
        return 1;
    }

    pairing_init_set_buf(pairing, param, count);

    std::printf("blocks,dec_keygen_ms,sign_keygen_ms,encryption_ms,signature_ms,decryption_ms,verification_ms\n");

    for (int bi = 0; bi < block_size_count; bi++) {
        int blocks = block_sizes[bi];
        PublicContext ctx;
        LsssPolicy policy;
        AuditorKey ask;

        init_public_context(&ctx, pairing, params);
        init_policy(&policy, pairing, params.l_rows, params.n_attrs, params.n_attrs);

        FileCiphertext* cts = static_cast<FileCiphertext*>(std::calloc(static_cast<size_t>(params.file_count), sizeof(FileCiphertext)));
        SigningKey* sks = static_cast<SigningKey*>(std::calloc(static_cast<size_t>(params.file_count), sizeof(SigningKey)));
        BlockSignature** sigs = static_cast<BlockSignature**>(std::calloc(static_cast<size_t>(params.file_count), sizeof(BlockSignature*)));
        die_alloc(cts);
        die_alloc(sks);
        die_alloc(sigs);
        for (int j = 0; j < params.file_count; j++) {
            sigs[j] = static_cast<BlockSignature*>(std::calloc(static_cast<size_t>(blocks), sizeof(BlockSignature)));
            die_alloc(sigs[j]);
        }

        double dec_keygen_ms = time_dec_keygen(&ctx, &ask);
        double encryption_ms = time_encrypt(&ctx, &policy, cts);
        double sign_keygen_ms = time_sign_keygen(&ctx, &policy, cts, sks);
        double signature_ms = time_signature(&ctx, &policy, cts, sks, sigs, blocks);
        double decryption_ms = time_decryption(&ctx, &policy, &ask, cts);
        double verification_ms = time_verification(&ctx, cts, sigs, blocks);

        std::printf("%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            blocks,
            dec_keygen_ms,
            sign_keygen_ms,
            encryption_ms,
            signature_ms,
            decryption_ms,
            verification_ms);

        for (int j = 0; j < params.file_count; j++) {
            clear_block_signatures(sigs[j], blocks);
            std::free(sigs[j]);
            clear_signing_key(&sks[j], &policy, params.n_attrs);
            clear_file_ciphertext(&cts[j], &policy);
        }
        clear_auditor_key(&ask);
        clear_policy(&policy);
        clear_public_context(&ctx);
        std::free(sigs);
        std::free(sks);
        std::free(cts);
    }

    pairing_clear(pairing);
    return 0;
}
