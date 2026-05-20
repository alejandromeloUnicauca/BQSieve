#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "polynomial.h"
#include "structsqs.h"


// Declaración de shanksTonelli usada en sieve.c para calcular raíces mod p
int shanksTonelli(mpz_t n, mpz_t p, mpz_t r1, mpz_t r2);

/*--------------------------------------------------------------------
 * SIQS polynomial generation: multiple b per each a
 *
 * a = q_0 * q_1 * ... * q_{s-1}  (product of s small primes from the factor base)
 * For each factor q_j, compute:
 *   B_j = (a/q_j) * ((a/q_j)^{-1} mod q_j) * (sqrt(N) mod q_j)
 * First b = sum(B_j), subsequent b's via Gray code: b_{i} = b_{i-1} ± 2*B_j
 * Total: 2^(s-1) polynomials for each 'a'
 *
 * The quadratic: Q(x) = a*x^2 + 2*b*x + c  where c = (b^2 - N)/a
 * Relation: (a*x + b)^2 = a * Q(x) + N  =>  (a*x+b)^2 ≡ a*Q(x) (mod N)
 *
 * For the square root phase we store lhs = a*x+b and Qval = Q(x).
 * Since lhs^2 = a*Q(x) + N, we need roota = sqrt(a) so that
 * lhs^2 = roota^2 * Q(x) + N => lhs^2 ≡ roota^2 * Q(x) (mod N).
 * But a is not a perfect square. Instead we note:
 *   For each relation: (a*x+b)^2 ≡ a*Q(x) (mod N)
 * The product over a subset S: prod((a_i*x_i+b_i)^2) ≡ prod(a_i*Q(x_i)) (mod N)
 * We need prod(a_i*Q(x_i)) to be a perfect square.
 * The matrix tracks parities of exponents in Q(x_i), plus the factors of a_i.
 * Since a_i's factors are primes from the base, their exponents are tracked in the matrix.
 *
 * Implementation: the trial division code in factoring.c divides Q(x) by primes
 * in the base. The factors of 'a' also divide Q(x) (because Q(x) = (a*x+b)^2/a - N/a
 * and b ≡ sqrt(N) mod q_j for each factor q_j of a, so q_j | Q(x) for x=0, etc.)
 * Actually, the factors of 'a' divide a*Q(x) = (a*x+b)^2 - N. Whether they divide
 * Q(x) = (a*x+b)^2/a - N/a... they do because b^2 ≡ N (mod a), so c = (b^2-N)/a is integer.
 *
 * Key insight: since a is composed of primes from the factor base, and trial division
 * on Q(x) will find those primes (they divide Q(x) at certain x), the exponent tracking
 * handles automatically. But we also need to track 'a' in the relation written to
 * polinomio.txt. We set roota = 1 and write Q_file = Q(x) (not a*Q(x)).
 * Then in mulPoli: prod(roota^2 * Q_i) = prod(Q_i). For the relation:
 *   lhs^2 = a*Q(x) + N => lhs^2 - N = a*Q(x)
 * This means lhs^2 ≢ Q(x) (mod N), but rather lhs^2 ≡ a*Q(x) (mod N).
 *
 * FIX: we need to multiply Q(x) by 'a' in the file, OR include 'a' as part of the
 * square root. The simplest: write Q_file = a*Q(x) = (a*x+b)^2 - N, set roota = 1.
 * Then lhs^2 = Q_file + N, so lhs^2 ≡ Q_file (mod N). Perfect.
 * Trial division then operates on Q_file/a = Q(x), but the exponent vector
 * must also include the factors of 'a'. We'll add them explicitly.
 *--------------------------------------------------------------------*/

/**
 * @brief Determina cuántos factores tendrá 'a' según el tamaño en bits de target_a.
 */
static unsigned int choose_num_factors(unsigned int a_bits) {
    if (a_bits <= 30)  return 2;
    if (a_bits <= 50)  return 3;
    if (a_bits <= 70)  return 4;
    if (a_bits <= 100) return 5;
    if (a_bits <= 130) return 6;
    if (a_bits <= 160) return 7;
    if (a_bits <= 200) return 8;
    return 9;
}

/**
 * @brief Selecciona los factores primos de 'a' desde la base de primos.
 */
static void select_a_factors(qs_struct *qs_data, siqs_poly_state *st) {
    unsigned int a_bits = mpz_sizeinbase(st->target_a, 2);
    unsigned int s = choose_num_factors(a_bits);
    if (s > MAX_SIQS_FACTORS) s = MAX_SIQS_FACTORS;
    if (s < 2) s = 2;
    st->num_factors = s;

    long base_len = qs_data->base.length;

    /* Rango de selección: usar primos del 20%-80% de la base */
    long pool_start = base_len / 5;
    long pool_end = (base_len * 4) / 5;
    if (pool_start < 3) pool_start = 3;
    if (pool_end <= pool_start + (long)s) pool_end = base_len - 1;
    long pool_size = pool_end - pool_start;
    if (pool_size < (long)s) {
        pool_start = 3;
        pool_end = base_len - 1;
        pool_size = pool_end - pool_start;
    }

    /* Seleccionar s-1 factores aleatoriamente (sin repetir) */
    mpz_t a_partial;
    mpz_init(a_partial);
    mpz_set_ui(a_partial, 1);

    for (unsigned int i = 0; i < s - 1; i++) {
        int unique;
        long idx;
        do {
            unique = 1;
            idx = pool_start + (random() % pool_size);
            for (unsigned int j = 0; j < i; j++) {
                if (st->factor_fb_idx[j] == (unsigned long)idx) {
                    unique = 0;
                    break;
                }
            }
        } while (!unique);

        st->factor_fb_idx[i] = (unsigned long)idx;
        mpz_set(st->factors[i], qs_data->base.primes[idx].value);
        mpz_mul(a_partial, a_partial, st->factors[i]);
    }

    /* Último factor: elegir el que acerque a_partial * q_{s-1} a target_a */
    mpz_t ideal, best_dist, dist;
    mpz_inits(ideal, best_dist, dist, NULL);
    mpz_fdiv_q(ideal, st->target_a, a_partial);

    long best_idx = -1;
    mpz_set_ui(best_dist, 0);
    mpz_setbit(best_dist, 128);

    for (long k = pool_start; k < pool_end; k++) {
        int dup = 0;
        for (unsigned int j = 0; j < s - 1; j++) {
            if (st->factor_fb_idx[j] == (unsigned long)k) { dup = 1; break; }
        }
        if (dup) continue;

        mpz_sub(dist, qs_data->base.primes[k].value, ideal);
        mpz_abs(dist, dist);
        if (mpz_cmp(dist, best_dist) < 0) {
            mpz_set(best_dist, dist);
            best_idx = k;
        }
    }
    mpz_clears(ideal, best_dist, dist, NULL);

    if (best_idx < 0) best_idx = pool_start;

    st->factor_fb_idx[s - 1] = (unsigned long)best_idx;
    mpz_set(st->factors[s - 1], qs_data->base.primes[best_idx].value);

    mpz_clear(a_partial);
}

/**
 * @brief Genera un nuevo valor de 'a' y calcula todos los B_j auxiliares.
 */
static void build_new_a(qs_struct *qs_data) {
    siqs_poly_state *st = &qs_data->siqs_state;

    select_a_factors(qs_data, st);
    unsigned int s = st->num_factors;

    /* a = producto de todos los factores */
    mpz_set_ui(qs_data->poly.a, 1);
    for (unsigned int i = 0; i < s; i++) {
        mpz_mul(qs_data->poly.a, qs_data->poly.a, st->factors[i]);
    }

    /* Para cada factor q_j, calcular B_j */
    mpz_t a_div_qj, gamma, sqrt_n_mod_qj, r2_tmp;
    mpz_inits(a_div_qj, gamma, sqrt_n_mod_qj, r2_tmp, NULL);

    mpz_set_ui(qs_data->poly.b, 0);

    for (unsigned int j = 0; j < s; j++) {
        mpz_t qj;
        mpz_init_set(qj, st->factors[j]);

        mpz_divexact(a_div_qj, qs_data->poly.a, qj);

        shanksTonelli(qs_data->n, qj, sqrt_n_mod_qj, r2_tmp);

        mpz_t inv_a_div_qj;
        mpz_init(inv_a_div_qj);
        mpz_mod(gamma, a_div_qj, qj);
        if (mpz_invert(inv_a_div_qj, gamma, qj) == 0) {
            mpz_set(sqrt_n_mod_qj, r2_tmp);
            mpz_mod(gamma, a_div_qj, qj);
            mpz_invert(inv_a_div_qj, gamma, qj);
        }
        mpz_mul(gamma, inv_a_div_qj, sqrt_n_mod_qj);
        mpz_mod(gamma, gamma, qj);

        /* Si gamma > q_j/2, gamma = q_j - gamma */
        mpz_t half_qj;
        mpz_init(half_qj);
        mpz_fdiv_q_ui(half_qj, qj, 2);
        if (mpz_cmp(gamma, half_qj) > 0) {
            mpz_sub(gamma, qj, gamma);
        }
        mpz_clear(half_qj);

        /* B_j = a_div_qj * gamma */
        mpz_mul(st->Bvals[j], a_div_qj, gamma);

        mpz_add(qs_data->poly.b, qs_data->poly.b, st->Bvals[j]);

        mpz_clear(inv_a_div_qj);
        mpz_clear(qj);
    }

    mpz_clears(a_div_qj, gamma, sqrt_n_mod_qj, r2_tmp, NULL);

    /* c = (b^2 - N) / a */
    mpz_mul(qs_data->poly.c, qs_data->poly.b, qs_data->poly.b);
    mpz_sub(qs_data->poly.c, qs_data->poly.c, qs_data->n);
    mpz_divexact(qs_data->poly.c, qs_data->poly.c, qs_data->poly.a);

    /* roota = 1 para SIQS (la relación se maneja con a*Q(x)) */
    mpz_set_ui(qs_data->roota, 1);

    /* Inicializar iteración de polinomios derivados */
    st->poly_index = 0;
    st->num_derived = 1UL << (s - 1);
    st->sieve_new_a = 1; /* la criba debe recomputar raíces y deltas */

    /* Doblar los B_j para futuras iteraciones */
    for (unsigned int j = 0; j < s; j++) {
        mpz_mul_ui(st->Bvals[j], st->Bvals[j], 2);
    }
}

/**
 * @brief Genera el siguiente polinomio derivado (nuevo b para el mismo a).
 * @return 1 si se generó, 0 si se agotaron los b's para este a
 */
static int next_siqs_b(qs_struct *qs_data) {
    siqs_poly_state *st = &qs_data->siqs_state;

    if (st->poly_index >= st->num_derived)
        return 0;

    if (st->poly_index == 0) {
        st->poly_index++;
        return 1;
    }

    /* Gray code: bit que cambia */
    unsigned long i = st->poly_index;
    unsigned int j = 0;
    while ((i & (1UL << j)) == 0)
        j++;

    /* Sumar o restar 2*B_j. Señalizar a la criba el factor y signo
     * para que actualice las raíces incrementalmente (sin recomputar). */
    if (i & (1UL << (j + 1))) {
        mpz_add(qs_data->poly.b, qs_data->poly.b, st->Bvals[j]);
        st->sieve_flip_sign = +1;
    } else {
        mpz_sub(qs_data->poly.b, qs_data->poly.b, st->Bvals[j]);
        st->sieve_flip_sign = -1;
    }
    st->sieve_flip_j = j;
    st->sieve_new_a = 0;

    /* Recalcular c = (b² - N) / a */
    mpz_mul(qs_data->poly.c, qs_data->poly.b, qs_data->poly.b);
    mpz_sub(qs_data->poly.c, qs_data->poly.c, qs_data->n);
    mpz_divexact(qs_data->poly.c, qs_data->poly.c, qs_data->poly.a);

    st->poly_index++;
    return 1;
}

/**
 * @brief Punto de entrada para generar el siguiente polinomio SIQS.
 */
int generate_mpqs_poly(qs_struct *qs_data) {
    siqs_poly_state *st = &qs_data->siqs_state;

    if (!st->initialized) {
        for (unsigned int i = 0; i < MAX_SIQS_FACTORS; i++) {
            mpz_init(st->factors[i]);
            mpz_init(st->Bvals[i]);
        }
        mpz_init(st->target_a);

        /* target_a = sqrt(2*N) / sieve_size */
        mpz_t root2n;
        mpz_init(root2n);
        mpz_mul_ui(root2n, qs_data->n, 2);
        mpz_sqrt(root2n, root2n);
        unsigned long ss = qs_data->sieve_params.sieve_size;
        if (ss == 0) ss = 65536;
        mpz_fdiv_q_ui(st->target_a, root2n, ss);
        mpz_clear(root2n);

        srandom((unsigned int)mpz_fdiv_ui(qs_data->n, 1000000007) ^ 0xdeadbeef);

        st->initialized = 1;
        st->poly_index = 0;
        st->num_derived = 0;
    }

    /* Intentar siguiente b derivado, si no hay más, nuevo a */
    if (!next_siqs_b(qs_data)) {
        build_new_a(qs_data);
        next_siqs_b(qs_data);
    }

    return 1;
}

/**
 * @brief Evalúa Q(x) = a*x² + 2*b*x + c para el polinomio actual.
 *
 * Este es el valor que los primos de la base dividen. Para la relación
 * (a*x+b)² ≡ a*Q(x) (mod N), el archivo polinomio.txt guardará Q(x)
 * y el factor 'a' se maneja vía la columna de exponentes de los factores de a.
 */
int eval_mpqs_Qx(qs_struct *qs_data, mpz_t x, mpz_t result){
    mpz_t ax, axx, two_bx, tmp;
    mpz_inits(ax, axx, two_bx, tmp, NULL);

    mpz_mul(ax, qs_data->poly.a, x);
    mpz_mul(axx, ax, x);
    mpz_mul(two_bx, qs_data->poly.b, x);
    mpz_mul_ui(two_bx, two_bx, 2);
    mpz_add(tmp, axx, two_bx);
    mpz_add(result, tmp, qs_data->poly.c);

    mpz_clears(ax, axx, two_bx, tmp, NULL);
    return 1;
}

