#include <gmp.h>
#include <ecm.h>
#include <stdio.h>
int main() {
    mpz_t N, f;
    mpz_inits(N, f, NULL);
    mpz_set_str(N, "2821457049646560553430084268822185147537933587", 10);
    ecm_params params;
    ecm_init(params);
    params->sigma_is_A = -1; /* use sigma */
    mpz_set_ui(params->sigma, 42);
    int ret = ecm_factor(f, N, 1e6, params);
    if (ret > 0) {
        gmp_printf("Factor: %Zd\n", f);
        mpz_t q; mpz_init(q);
        mpz_tdiv_q(q, N, f);
        gmp_printf("Cofactor: %Zd\n", q);
        mpz_clear(q);
    } else {
        printf("ECM didn't find factor with B1=1e6\n");
    }
    ecm_clear(params);
    mpz_clears(N, f, NULL);
    return 0;
}
