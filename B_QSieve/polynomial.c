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
    // target_a = sqrt(2*N) / sieve_size  (valor óptimo de 'a')
    // Como a = roota², queremos roota ≈ target_a^(1/2) = (sqrt(2*N) / sieve_size)^(1/2)
    if (mpz_cmp_ui(qs_data->roota, 0) == 0) {
        mpz_t root2n, target_a;
        mpz_inits(root2n, target_a, NULL);

        // root2n = floor(sqrt(2*N))
        mpz_mul_ui(root2n, qs_data->n, 2);
        mpz_sqrt(root2n, root2n);

        // target_a = root2n / sieve_size
        unsigned long sieve_size = qs_data->sieve_params.sieve_size;
        if (sieve_size == 0) sieve_size = 65536; // fallback
        mpz_fdiv_q_ui(target_a, root2n, sieve_size);

        // roota = floor(sqrt(target_a))
        mpz_sqrt(qs_data->roota, target_a);
        if (mpz_cmp_ui(qs_data->roota, 3) < 0)
            mpz_set_ui(qs_data->roota, 3);
        // make odd if even
        if (mpz_even_p(qs_data->roota))
            mpz_add_ui(qs_data->roota, qs_data->roota, 1);

        mpz_clears(root2n, target_a, NULL);
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

