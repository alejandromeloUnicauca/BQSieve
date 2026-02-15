#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "polynomial.h"
#include "structsqs.h"


// Declaración de shanksTonelli usada en sieve.c para calcular raíces mod p
int shanksTonelli(mpz_t n, mpz_t p, mpz_t r1, mpz_t r2);

int generate_mpqs_poly(qs_struct *qs_data){
    // Usa qs_data->roota persistente entre llamadas.
    // Primera llamada (roota==0): calcula valor inicial.
    // Cada llamada: avanza roota a next_prime con legendre==1.

    mpz_t a_mp, b_mp, c_mp, inv_mp, two_b_mod;
    mpz_inits(a_mp, b_mp, c_mp, inv_mp, two_b_mod, NULL);
    mpz_t qtmp, tmp2;
    mpz_inits(qtmp, tmp2, NULL);

    // Si roota == 0, calcular valor inicial
    if (mpz_cmp_ui(qs_data->roota, 0) == 0) {
        // 1) calcular bound = int(5 * (log10(n))^2)
        mpfr_t nn, log10n, tmpfr;
        mpfr_inits(nn, log10n, tmpfr, NULL);
        mpfr_set_z(nn, qs_data->n, MPFR_RNDN);
        mpfr_log(tmpfr, nn, MPFR_RNDZ);
        mpfr_set_str(log10n, "2.302585092994046", 10, MPFR_RNDZ);
        mpfr_div(log10n, tmpfr, log10n, MPFR_RNDZ);
        mpfr_mul(tmpfr, log10n, log10n, MPFR_RNDZ);
        mpfr_mul_ui(tmpfr, tmpfr, 5, MPFR_RNDZ);
        unsigned long bound = mpfr_get_ui(tmpfr, MPFR_RNDZ);
        mpfr_clears(nn, log10n, tmpfr, NULL);

        // construir factorbase temporal contando primos <= bound con legendre==1
        unsigned long fb_count = 0;
        for (unsigned long p = 2; p <= bound; p++) {
            int is_prime = 1;
            if (p < 2) is_prime = 0;
            for (unsigned long d = 2; d * d <= p && is_prime; d++) {
                if (p % d == 0) { is_prime = 0; }
            }
            if (!is_prime) continue;
            mpz_t p_mp;
            mpz_init_set_ui(p_mp, p);
            int leg = mpz_legendre(qs_data->n, p_mp);
            mpz_clear(p_mp);
            if (leg == 1 || p == 2) {
                fb_count++;
            }
        }

        unsigned long xmax = fb_count * 60 * 4;
        if (xmax == 0) xmax = 1;

        // calcular root2n = floor(sqrt(2*n))
        mpz_t root2n, tmpz;
        mpz_inits(root2n, tmpz, NULL);
        mpz_mul_ui(root2n, qs_data->n, 2);
        mpz_sqrt(root2n, root2n);

        // roota = isqrt(root2n // xmax)
        mpz_fdiv_q_ui(tmpz, root2n, xmax);
        mpz_sqrt(qs_data->roota, tmpz);
        if (mpz_cmp_ui(qs_data->roota, 3) < 0)
            mpz_set_ui(qs_data->roota, 3);
        // make odd if even
        if (mpz_even_p(qs_data->roota))
            mpz_add_ui(qs_data->roota, qs_data->roota, 1);

        mpz_clears(root2n, tmpz, NULL);
    }

    // Avanzar roota al siguiente primo con legendre(n, roota) == 1
    // (equivalente al loop Python: while 1: roota=next_prime(roota); if legendre==1: break)
    int found = 0;
    while (!found) {
        mpz_nextprime(qs_data->roota, qs_data->roota);
        if (mpz_legendre(qs_data->n, qs_data->roota) == 1) {
            found = 1;
        }
    }

    // a = roota^2
    mpz_mul(a_mp, qs_data->roota, qs_data->roota);

    // compute b = modular sqrt n mod roota
    mpz_t r1, r2;
    mpz_inits(r1, r2, NULL);
    shanksTonelli(qs_data->n, qs_data->roota, r1, r2);
    mpz_set(b_mp, r1);

    // compute inv of 2*b mod roota
    mpz_mul_ui(qtmp, b_mp, 2);
    mpz_mod(two_b_mod, qtmp, qs_data->roota);
    if (mpz_invert(inv_mp, two_b_mod, qs_data->roota) == 0) {
        // fallback: try r2 as b
        mpz_set(b_mp, r2);
        mpz_mul_ui(qtmp, b_mp, 2);
        mpz_mod(two_b_mod, qtmp, qs_data->roota);
        if (mpz_invert(inv_mp, two_b_mod, qs_data->roota) == 0) {
            // give up: use simple roota=3 fallback
            mpz_set_ui(qs_data->roota, 3);
            mpz_mul(a_mp, qs_data->roota, qs_data->roota);
            mpz_sqrt(qtmp, qs_data->n);
            mpz_mod(b_mp, qtmp, a_mp);
        }
    }

    // b = (b - (b*b - N) * inv) mod a
    mpz_mul(qtmp, b_mp, b_mp);
    mpz_sub(qtmp, qtmp, qs_data->n);
    mpz_mul(tmp2, qtmp, inv_mp);
    mpz_sub(qtmp, b_mp, tmp2);
    mpz_mod(b_mp, qtmp, a_mp);

    // c = (b^2 - N) / a
    mpz_mul(qtmp, b_mp, b_mp);
    mpz_sub(qtmp, qtmp, qs_data->n);
    if (mpz_divisible_p(qtmp, a_mp)) {
        mpz_divexact(c_mp, qtmp, a_mp);
    } else {
        mpz_fdiv_q(c_mp, qtmp, a_mp);
    }

    // assign to qs_data->poly
    mpz_set(qs_data->poly.a, a_mp);
    mpz_set(qs_data->poly.b, b_mp);
    mpz_set(qs_data->poly.c, c_mp);

    // cleanup
    mpz_clears(a_mp, b_mp, c_mp, inv_mp, two_b_mod, qtmp, tmp2, r1, r2, NULL);

    return 1;
}

int eval_mpqs_Qx(qs_struct *qs_data, mpz_t x, mpz_t result){
    // result = a*x*x + 2*b*x + c  (same definition used in the Python POC)
    mpz_t ax, axx, two_bx, tmp;
    mpz_inits(ax, axx, two_bx, tmp, NULL);

    // ax = a * x
    mpz_mul(ax, qs_data->poly.a, x);
    // axx = ax * x = a * x * x
    mpz_mul(axx, ax, x);
    // two_bx = 2 * b * x
    mpz_mul(two_bx, qs_data->poly.b, x);
    mpz_mul_ui(two_bx, two_bx, 2);
    // tmp = axx + two_bx
    mpz_add(tmp, axx, two_bx);
    // result = tmp + c
    mpz_add(result, tmp, qs_data->poly.c);

    mpz_clears(ax, axx, two_bx, tmp, NULL);
    return 1;
}

