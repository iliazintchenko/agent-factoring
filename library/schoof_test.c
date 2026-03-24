#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

/* Simple polynomial over Z/NZ */
typedef struct { mpz_t *c; int deg, alloc; } poly_t;

void p_init(poly_t *p, int n) {
    p->c = malloc(n * sizeof(mpz_t));
    for (int i = 0; i < n; i++) mpz_init(p->c[i]);
    p->deg = -1; p->alloc = n;
}
void p_free(poly_t *p) {
    for (int i = 0; i < p->alloc; i++) mpz_clear(p->c[i]);
    free(p->c);
}

/* Multiply a*b mod (mod_poly, N), store in res. Returns 1 if factor found. */
int p_mulmod(poly_t *res, const poly_t *a, const poly_t *b,
             const poly_t *m, const mpz_t N, mpz_t factor) {
    if (a->deg < 0 || b->deg < 0) { res->deg = -1; return 0; }
    int pd = a->deg + b->deg;
    int md = m->deg;
    
    /* Temp storage */
    int sz = pd + 1;
    mpz_t *t = malloc(sz * sizeof(mpz_t));
    for (int i = 0; i < sz; i++) { mpz_init(t[i]); mpz_set_ui(t[i], 0); }
    
    /* Multiply */
    for (int i = 0; i <= a->deg; i++)
        for (int j = 0; j <= b->deg; j++) {
            mpz_addmul(t[i+j], a->c[i], b->c[j]);
            mpz_mod(t[i+j], t[i+j], N);
        }
    
    /* Reduce mod m */
    mpz_t lci, q;
    mpz_init(lci); mpz_init(q);
    if (!mpz_invert(lci, m->c[md], N)) {
        mpz_gcd(factor, m->c[md], N);
        mpz_clear(lci); mpz_clear(q);
        for (int i = 0; i < sz; i++) mpz_clear(t[i]); free(t);
        return 1;
    }
    for (int i = pd; i >= md; i--) {
        if (mpz_sgn(t[i]) == 0) continue;
        mpz_mul(q, t[i], lci); mpz_mod(q, q, N);
        for (int j = 0; j <= md; j++) {
            mpz_submul(t[i - md + j], q, m->c[j]);
            mpz_mod(t[i-md+j], t[i-md+j], N);
        }
    }
    
    /* Copy to res (degree < md) */
    res->deg = md - 1;
    while (res->deg >= 0 && mpz_sgn(t[res->deg]) == 0) res->deg--;
    for (int i = 0; i <= res->deg && i < res->alloc; i++)
        mpz_set(res->c[i], t[i]);
    
    mpz_clear(lci); mpz_clear(q);
    for (int i = 0; i < sz; i++) mpz_clear(t[i]); free(t);
    return 0;
}

/* x^exp mod (m, N) */
int p_powx(poly_t *res, const mpz_t exp, const poly_t *m, const mpz_t N, mpz_t factor) {
    /* res = 1 */
    res->deg = 0; mpz_set_ui(res->c[0], 1);
    
    /* base = x */
    poly_t base; p_init(&base, m->alloc);
    base.deg = 1; mpz_set_ui(base.c[0], 0); mpz_set_ui(base.c[1], 1);
    
    poly_t tmp; p_init(&tmp, m->alloc);
    
    size_t bits = mpz_sizeinbase(exp, 2);
    for (size_t i = 0; i < bits; i++) {
        if (mpz_tstbit(exp, i)) {
            if (p_mulmod(&tmp, res, &base, m, N, factor)) {
                p_free(&base); p_free(&tmp); return 1;
            }
            for (int k = 0; k <= tmp.deg && k < res->alloc; k++)
                mpz_set(res->c[k], tmp.c[k]);
            for (int k = tmp.deg + 1; k < res->alloc; k++) mpz_set_ui(res->c[k], 0);
            res->deg = tmp.deg;
        }
        if (i + 1 < bits) {
            if (p_mulmod(&tmp, &base, &base, m, N, factor)) {
                p_free(&base); p_free(&tmp); return 1;
            }
            for (int k = 0; k <= tmp.deg && k < base.alloc; k++)
                mpz_set(base.c[k], tmp.c[k]);
            for (int k = tmp.deg + 1; k < base.alloc; k++) mpz_set_ui(base.c[k], 0);
            base.deg = tmp.deg;
        }
    }
    p_free(&base); p_free(&tmp);
    return 0;
}

int main(int argc, char **argv) {
    mpz_t N, factor;
    mpz_init_set_ui(N, 143); /* 11 * 13 */
    mpz_init(factor);
    
    gmp_randstate_t rng;
    gmp_randinit_mt(rng);
    gmp_randseed_ui(rng, 42);
    
    int found = 0;
    for (int curve = 0; curve < 100 && !found; curve++) {
        mpz_t a, b;
        mpz_init(a); mpz_init(b);
        mpz_urandomm(a, rng, N);
        mpz_urandomm(b, rng, N);
        
        /* mod = x^3 + ax + b */
        poly_t mod; p_init(&mod, 4);
        mpz_set(mod.c[0], b); mpz_set(mod.c[1], a);
        mpz_set_ui(mod.c[2], 0); mpz_set_ui(mod.c[3], 1);
        mod.deg = 3;
        
        /* Compute x^N mod (x^3+ax+b, N) */
        poly_t res; p_init(&res, 4);
        fprintf(stderr, "curve %d: a=%s b=%s\n", curve, mpz_get_str(NULL,10,a), mpz_get_str(NULL,10,b));
        
        int ret = p_powx(&res, N, &mod, N, factor);
        if (ret) {
            if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
                gmp_printf("FACTOR: %Zd (curve %d)\n", factor, curve);
                found = 1;
            }
        } else {
            /* Check gcd(x^N - x, x^3+ax+b) mod N */
            if (res.deg >= 1) mpz_sub_ui(res.c[1], res.c[1], 1);
            else { res.deg = 1; mpz_sub_ui(res.c[1], N, 1); }
            mpz_mod(res.c[1], res.c[1], N);
            
            /* Check coefficients for factors */
            for (int i = 0; i <= res.deg; i++) {
                mpz_t g; mpz_init(g);
                mpz_gcd(g, res.c[i], N);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                    gmp_printf("FACTOR: %Zd (curve %d, coef %d)\n", g, curve, i);
                    found = 1;
                }
                mpz_clear(g);
                if (found) break;
            }
        }
        
        p_free(&mod); p_free(&res);
        mpz_clear(a); mpz_clear(b);
    }
    
    if (!found) fprintf(stderr, "No factor found\n");
    
    gmp_randclear(rng);
    mpz_clear(N); mpz_clear(factor);
    return found ? 0 : 1;
}
