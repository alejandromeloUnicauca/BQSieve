#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include <math.h>
#include <omp.h>
#include "structsqs.h"

FILE *g_polinomio_fp = NULL;

void polinomio_open(void) {
    if (g_polinomio_fp == NULL) {
        g_polinomio_fp = fopen("polinomio.txt", "w");
        if (g_polinomio_fp == NULL) {
            perror("fopen polinomio.txt");
            exit(EXIT_FAILURE);
        }
    }
}

void polinomio_close(void) {
    if (g_polinomio_fp != NULL) {
        fclose(g_polinomio_fp);
        g_polinomio_fp = NULL;
    }
}

static double g_t_trialDiv_total = 0;
static double g_t_calloc         = 0;
static double g_t_recip_loop     = 0;
static double g_t_gmp_div        = 0;
static double g_t_prime_test     = 0;
static double g_t_combine        = 0;
static double g_t_pollard        = 0;
static double g_t_insert_matrix  = 0;
static double g_t_writeFullRel   = 0;
/* contadores específicos del path blockDivisionV2 */
static double g_t_block_gcd      = 0;
static double g_t_block_divexact = 0;
static double g_t_block_extract  = 0;
static unsigned long g_n_calls       = 0;
static unsigned long g_n_full        = 0;
static unsigned long g_n_1lp_try     = 0;
static unsigned long g_n_2lp_try     = 0;
static unsigned long g_n_pollard_run = 0;
static unsigned long g_n_writes      = 0;
double get_writeFullRel_time(void) { return g_t_writeFullRel; }
void   add_writeFullRel_time(double dt) {
    #pragma omp atomic update
    g_t_writeFullRel += dt;
    #pragma omp atomic update
    g_n_writes++;
}

void print_factoring_stats(double total_wall) {
    if (total_wall <= 0) total_wall = 1e-9;
    double sum_trial = g_t_calloc + g_t_recip_loop + g_t_gmp_div;
    double sum_block = g_t_block_gcd + g_t_block_divexact + g_t_block_extract;
    double sum_classify = g_t_prime_test + g_t_combine + g_t_pollard + g_t_insert_matrix;
    double other_in_td = g_t_trialDiv_total - sum_trial - sum_block - sum_classify;
    if (other_in_td < 0) other_in_td = 0;
    double outside_td = total_wall - g_t_trialDiv_total - g_t_writeFullRel;
    if (outside_td < 0) outside_td = 0;
    int is_blocks = (sum_block > 0);
    const char *fname = is_blocks ? "blockDivisionV2" : "trialDivisionRecip";
    printf("\n  [trial+combine breakdown]\n");
    printf("    candidatos procesados : %lu\n", g_n_calls);
    printf("    full relations        : %lu  (writes polinomio.txt: %lu)\n", g_n_full, g_n_writes);
    printf("    1LP candidatos        : %lu\n", g_n_1lp_try);
    printf("    2LP candidatos        : %lu  (Pollard rho corrido %lu veces)\n",
           g_n_2lp_try, g_n_pollard_run);
    printf("    %s total: %.3fs (%5.1f%%)\n", fname, g_t_trialDiv_total, 100.0*g_t_trialDiv_total/total_wall);
    if (is_blocks) {
        printf("      ┌─ mpz_gcd          : %.3fs (%5.1f%%)\n", g_t_block_gcd,      100.0*g_t_block_gcd/total_wall);
        printf("      ├─ mpz_divexact     : %.3fs (%5.1f%%)\n", g_t_block_divexact, 100.0*g_t_block_divexact/total_wall);
        printf("      ├─ extract exp_vec  : %.3fs (%5.1f%%)\n", g_t_block_extract,  100.0*g_t_block_extract/total_wall);
    } else {
        printf("      ┌─ calloc exp_vec   : %.3fs (%5.1f%%)\n", g_t_calloc,         100.0*g_t_calloc/total_wall);
        printf("      ├─ recip test loop  : %.3fs (%5.1f%%)\n", g_t_recip_loop,     100.0*g_t_recip_loop/total_wall);
        printf("      ├─ mpz_tdiv_q_ui    : %.3fs (%5.1f%%)\n", g_t_gmp_div,        100.0*g_t_gmp_div/total_wall);
    }
    printf("      ├─ mpz_probab_prime : %.3fs (%5.1f%%)\n", g_t_prime_test,    100.0*g_t_prime_test/total_wall);
    printf("      ├─ try_combine_part : %.3fs (%5.1f%%)\n", g_t_combine,       100.0*g_t_combine/total_wall);
    printf("      ├─ Pollard rho      : %.3fs (%5.1f%%)\n", g_t_pollard,       100.0*g_t_pollard/total_wall);
    printf("      ├─ insertarNumero(M): %.3fs (%5.1f%%)\n", g_t_insert_matrix, 100.0*g_t_insert_matrix/total_wall);
    printf("      └─ otros (mpz_init/clear, etc): %.3fs (%5.1f%%)\n",
           other_in_td, 100.0*other_in_td/total_wall);
    printf("    escritura polinomio.txt (full rels): %.3fs (%5.1f%%)\n", g_t_writeFullRel, 100.0*g_t_writeFullRel/total_wall);
    printf("    fuera de %s (loop overhead): %.3fs (%5.1f%%)\n", fname, outside_td, 100.0*outside_td/total_wall);
    fflush(stdout);
}

int trialDivisionRecip(mpz_t Qxi, qs_struct *qs_data, mpz_t Xi, unsigned long sieve_offset);
int blockDivisionV2(mpz_t Qxi, qs_struct *qs_data, mpz_t Xi);
static int classify_smooth_or_partial(qs_struct *qs_data, mpz_t Qxi, mpz_t Xi,
                                       int *exp_vec, int sign, mpz_t residual);
void insertarNumero(matrix * matriz, int posFila, int posColumna, int valor);
int factor_cofactor_pollard(mpz_t cofactor, unsigned long lp_bound,
                            unsigned long *f1, unsigned long *f2);
int try_combine_partial(qs_struct *qs_data, mpz_t Qxi, mpz_t Xi,
                        int *exp_vec, int sign,
                        unsigned long lp1, unsigned long lp2);

/**
 * @brief Añade los factores de 'a' al vector de exponentes en la matriz.
 *
 * En SIQS, a = q_0 * q_1 * ... * q_{s-1} donde cada q_j está en la base.
 * La relación es lhs² ≡ a*Q(x) (mod N), y el vector de exponentes de Q(x)
 * ya está en la matriz. Necesitamos sumar exponente 1 por cada factor de a.
 * factor_fb_idx[j] da el índice en la base → columna = idx+1 en la matriz.
 */
void add_a_factors_to_matrix(qs_struct *qs_data) {
    siqs_poly_state *st = &qs_data->siqs_state;
    for (unsigned int j = 0; j < st->num_factors; j++) {
        unsigned long fb_idx = st->factor_fb_idx[j];
        /* columna = fb_idx + 1 (columna 0 es signo) */
        insertarNumero(&qs_data->mat, qs_data->n_BSuaves, (int)fb_idx + 1, 1);
    }
}

void insertarNumero(matrix * matriz, int posFila, int posColumna, int valor){

	/*printf("Fila:%d Columna:%d Valor:%d\n",posFila,posColumna,valor);
	printf("Fila:%d Columna:%d Valor:%d\n",posFila,posColumna,matriz->data[posFila][posColumna]);
	fflush(stdout);*/
	
    // Verificar si la matriz y los índices son válidos antes de continuar
    if (matriz == NULL || matriz->data == NULL ||
        posFila < 0 || posFila >= matriz->n_rows ||	
        posColumna < 0 || posColumna >= matriz->n_cols) {
		fprintf(stderr,"Error insertarNumero: Fila=%d (max=%d) Columna=%d (max=%d) Valor=%d\n",
			posFila, matriz->n_rows, posColumna, matriz->n_cols, valor);
		fflush(stderr);
        exit(EXIT_FAILURE);
    }
	
	if(matriz->data != NULL){
		if(matriz->data != NULL && matriz->data[posFila][posColumna] == 0 && valor == 0){
			matriz->data[posFila][posColumna] = 0;
			return;
		}
		
		if(matriz->data[posFila][posColumna] == 1 && valor == 1){
			matriz->data[posFila][posColumna] = 0;
			return;
		}
		
		if(matriz->data[posFila][posColumna] == 0 && valor == 1){
			matriz->data[posFila][posColumna] = 1;
			return;
		}
		
		if(matriz->data[posFila][posColumna] == 1 && valor == 0){
			matriz->data[posFila][posColumna] = 1;
			return;
		}
	}
}

/*--------------------------------------------------------------------
 * blockDivisionV2 — Trial division agrupando primos por bloques.
 *
 * Por cada bloque k, computa gcd(Qxi, prod_factors[k]). Si gcd!=1,
 * divide y repite hasta que gcd vuelve a 1 — el bloque entero queda
 * procesado en O(K_iter × bloque) operaciones.
 *
 * Produce exp_vec[base.length] con los exponentes enteros y delega la
 * decisión full/1LP/2LP/discard a classify_smooth_or_partial.
 *
 * @return 1 = full, 2 = combined partial, 0 = saved/discarded.
 *--------------------------------------------------------------------*/
int blockDivisionV2(mpz_t Qxi, qs_struct *qs_data, mpz_t Xi)
{
    double _t_func0 = omp_get_wtime();
    #pragma omp atomic update
    g_n_calls++;

    int *exp_vec = (int *)calloc(qs_data->base.length, sizeof(int));
    mpz_t QxiTemp, gcd;
    mpz_inits(QxiTemp, gcd, NULL);
    mpz_abs(QxiTemp, Qxi);
    int sign = (mpz_sgn(Qxi) < 0) ? 1 : 0;

    long prime_offset = 0;
    for (long k = 0; k < (long)qs_data->blocks.length; k++) {
        long block_len = qs_data->blocks.block[k].length;

        double _tg0 = omp_get_wtime();
        mpz_gcd(gcd, QxiTemp, qs_data->blocks.block[k].prod_factors);
        #pragma omp atomic update
        g_t_block_gcd += omp_get_wtime() - _tg0;

        while (mpz_cmp_ui(gcd, 1) != 0) {
            double _td0 = omp_get_wtime();
            mpz_divexact(QxiTemp, QxiTemp, gcd);
            #pragma omp atomic update
            g_t_block_divexact += omp_get_wtime() - _td0;

            /* gcd es squarefree (prod_factors lo es), así que cada primo del
             * bloque que divide a gcd contribuye con exactamente +1 al
             * exponente en esta vuelta. */
            double _te0 = omp_get_wtime();
            for (long j = 0; j < block_len; j++) {
                uint32_t p = qs_data->base.primes[prime_offset + j].p;
                if (mpz_divisible_ui_p(gcd, p)) {
                    exp_vec[prime_offset + j]++;
                }
            }
            #pragma omp atomic update
            g_t_block_extract += omp_get_wtime() - _te0;

            double _tg1 = omp_get_wtime();
            mpz_gcd(gcd, QxiTemp, qs_data->blocks.block[k].prod_factors);
            #pragma omp atomic update
            g_t_block_gcd += omp_get_wtime() - _tg1;
        }

        prime_offset += block_len;
    }

    mpz_clear(gcd);

    int result = classify_smooth_or_partial(qs_data, Qxi, Xi, exp_vec, sign, QxiTemp);

    #pragma omp atomic update
    g_t_trialDiv_total += omp_get_wtime() - _t_func0;
    return result;
}

int factoringBlocks(qs_struct *qs_data, unsigned long endPos, unsigned long posXi, unsigned long xmax)
{
    (void)xmax;
    FILE *fp = g_polinomio_fp;

    for (unsigned long i = 0; i < endPos; i++) {
        int result = blockDivisionV2(qs_data->intervalo.Qxi[i], qs_data,
                                     qs_data->intervalo.Xi[posXi]);
        if (result == 1) {
            qs_data->n_BSuaves++;
            double _tw0 = omp_get_wtime();
            mpz_t lhs, Qfile;
            mpz_inits(lhs, Qfile, NULL);
            mpz_mul(lhs, qs_data->poly.a, qs_data->intervalo.Xi[posXi]);
            mpz_add(lhs, lhs, qs_data->poly.b);
            mpz_mul(Qfile, lhs, lhs);
            mpz_sub(Qfile, Qfile, qs_data->n);
            mpz_out_str(fp, 10, lhs);
            fprintf(fp, ";");
            mpz_out_str(fp, 10, Qfile);
            fprintf(fp, ";");
            mpz_out_str(fp, 10, qs_data->roota);
            fprintf(fp, "\n");
            mpz_clears(lhs, Qfile, NULL);
            add_writeFullRel_time(omp_get_wtime() - _tw0);
            if (qs_data->n_BSuaves == qs_data->mat.n_rows) return 0;
        } else if (result == 2) {
            qs_data->n_BSuaves++;
            if (qs_data->n_BSuaves == qs_data->mat.n_rows) return 0;
        }
        posXi++;
    }
    return 1;
}

/*--------------------------------------------------------------------
 * factor_cofactor_pollard — Factorizar un cofactor con Pollard-rho (GMP).
 *
 * Intenta dividir 'cofactor' en dos primos f1, f2 tales que
 * ambos < lp_bound. Usa la función mpz de GMP para Pollard-rho
 * con iteraciones limitadas.
 *
 * @param cofactor   Número a factorizar (compuesto, > 1)
 * @param lp_bound   Límite para large primes
 * @param f1, f2     Salida: los dos factores (si retorna 1)
 * @return 1 si se factorizó, 0 si no
 *--------------------------------------------------------------------*/
int factor_cofactor_pollard(mpz_t cofactor, unsigned long lp_bound,
                            unsigned long *f1, unsigned long *f2)
{
    /* Pollard-rho con función f(x) = x²+c mod n */
    mpz_t x, y, d, temp, n;
    mpz_inits(x, y, d, temp, n, NULL);
    mpz_set(n, cofactor);

    /* Probar varias semillas c */
    for (unsigned long c = 1; c <= 20; c++) {
        mpz_set_ui(x, 2);
        mpz_set_ui(y, 2);

        for (unsigned long iter = 0; iter < 10000; iter++) {
            /* x = x²+c mod n */
            mpz_mul(x, x, x);
            mpz_add_ui(x, x, c);
            mpz_mod(x, x, n);
            /* y = (y²+c)²+c mod n (paso doble) */
            mpz_mul(y, y, y);
            mpz_add_ui(y, y, c);
            mpz_mod(y, y, n);
            mpz_mul(y, y, y);
            mpz_add_ui(y, y, c);
            mpz_mod(y, y, n);
            /* d = gcd(|x-y|, n) */
            mpz_sub(temp, x, y);
            mpz_abs(temp, temp);
            mpz_gcd(d, temp, n);

            if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, n) < 0) {
                /* Factor encontrado */
                if (mpz_fits_ulong_p(d)) {
                    *f1 = mpz_get_ui(d);
                    mpz_divexact(temp, n, d);
                    if (mpz_fits_ulong_p(temp)) {
                        *f2 = mpz_get_ui(temp);
                        mpz_clears(x, y, d, temp, n, NULL);
                        return 1;
                    }
                }
                mpz_clears(x, y, d, temp, n, NULL);
                return 0;
            }
            if (mpz_cmp(d, n) == 0) break; /* ciclo, probar otra c */
        }
    }
    mpz_clears(x, y, d, temp, n, NULL);
    return 0;
}

/*--------------------------------------------------------------------
 * store_partial — Almacenar una relación parcial (1LP o 2LP)
 *--------------------------------------------------------------------*/
static void store_partial(qs_struct *qs_data, mpz_t Qxi, mpz_t Xi,
                          int *exp_vec, int sign,
                          unsigned long lp1, unsigned long lp2)
{
    if (qs_data->partials.n >= qs_data->partials.capacity) {
        unsigned long newcap = qs_data->partials.capacity == 0 ?
                               1024 : qs_data->partials.capacity * 2;
        qs_data->partials.entries = realloc(qs_data->partials.entries,
                                            newcap * sizeof(partial_entry));
        qs_data->partials.capacity = newcap;
    }
    unsigned long idx = qs_data->partials.n++;
    partial_entry *e = &qs_data->partials.entries[idx];
    e->large_prime = lp1;
    e->large_prime2 = lp2; /* 0 para 1LP, != 0 para 2LP */
    e->exponents = exp_vec; /* transferir ownership */
    e->sign = sign;
    mpz_init(e->lhs);
    mpz_mul(e->lhs, qs_data->poly.a, Xi);
    mpz_add(e->lhs, e->lhs, qs_data->poly.b);
    mpz_init_set(e->Qx, Qxi);
    mpz_init_set(e->roota, qs_data->roota);
    e->num_a_factors = qs_data->siqs_state.num_factors;
    for (unsigned int j = 0; j < e->num_a_factors; j++)
        e->a_factor_fb_idx[j] = qs_data->siqs_state.factor_fb_idx[j];
    mpz_init_set(e->a_value, qs_data->poly.a);
}

/*--------------------------------------------------------------------
 * combine_two_partials — Combinar dos parciales que comparten un primo grande.
 *
 * El primo compartido se cancela (exponente par → desaparece mod 2).
 * Si una o ambas son 2LP, los primos NO compartidos quedan como
 * new_lp1 y new_lp2 en la combinación. Si ambos no-compartidos
 * también se cancelan, la combinación es una full relation.
 *
 * Retorna: 1 = full relation directa,
 *          3 = nueva parcial (combinada pero con 1 o 2 LP residuales)
 *--------------------------------------------------------------------*/
static int combine_two_partials(qs_struct *qs_data,
                                int *exp_vec_a, int sign_a,
                                mpz_t Qxi_a, mpz_t Xi_a,
                                partial_entry *match,
                                unsigned long shared_lp,
                                unsigned long *remaining_lps,
                                int *n_remaining)
{
    (void)shared_lp; /* se cancela automáticamente en el XOR de exponentes */

    /* Recopilar los primos grandes NO compartidos de ambos lados.
     * lado A: exp_vec_a con primos lp1_a, lp2_a (ambos o uno provienen del caller)
     * lado B: match con lp1_b, lp2_b
     * El primo compartido (shared_lp) aparece en ambos lados → se cancela.
     * Los primos residuales son los que no se cancelan. */

    /* Escribir vector combinado (XOR) en la matriz */
    int combined_sign = (sign_a + match->sign) % 2;
    if (combined_sign)
        insertarNumero(&qs_data->mat, qs_data->n_BSuaves, 0, 1);
    for (long i = 0; i < qs_data->base.length; i++) {
        int combined_exp = (exp_vec_a[i] + match->exponents[i]) % 2;
        insertarNumero(&qs_data->mat, qs_data->n_BSuaves, i + 1, combined_exp);
    }
    add_a_factors_to_matrix(qs_data);
    for (unsigned int j = 0; j < match->num_a_factors; j++)
        insertarNumero(&qs_data->mat, qs_data->n_BSuaves,
                       (int)match->a_factor_fb_idx[j] + 1, 1);

    /* Escribir relación combinada en polinomio.txt */
    FILE *fp = g_polinomio_fp;
    if (fp) {
        mpz_t combined_lhs, combined_Q, my_lhs, aQ1, aQ2;
        mpz_inits(combined_lhs, combined_Q, my_lhs, aQ1, aQ2, NULL);
        mpz_mul(my_lhs, qs_data->poly.a, Xi_a);
        mpz_add(my_lhs, my_lhs, qs_data->poly.b);
        mpz_mul(combined_lhs, my_lhs, match->lhs);
        mpz_mul(aQ1, qs_data->poly.a, Qxi_a);
        mpz_mul(aQ2, match->a_value, match->Qx);
        mpz_mul(combined_Q, aQ1, aQ2);
        mpz_out_str(fp, 10, combined_lhs);
        fprintf(fp, ";");
        mpz_out_str(fp, 10, combined_Q);
        fprintf(fp, ";");
        mpz_out_str(fp, 10, qs_data->roota);
        fprintf(fp, ",");
        mpz_out_str(fp, 10, match->roota);
        fprintf(fp, "\n");
        mpz_clears(combined_lhs, combined_Q, my_lhs, aQ1, aQ2, NULL);
    }

    *n_remaining = 0;
    (void)remaining_lps;
    /* Nota: los primos no compartidos se tratan como factores del cofactor
     * combinado. Para que la relación sea full, necesitamos que TODOS los
     * primos grandes se cancelen. Para 1LP + 1LP con mismo LP → full.
     * Para 2LP + 1LP o 2LP + 2LP con un solo LP compartido → hay residuo.
     * En esta primera versión, solo combinamos 1LP con 1LP (ambas comparten
     * el único LP) y 2LP con 2LP que comparten ambos LPs.
     * Las combinaciones más complejas se dejan para futuras mejoras. */

    return 1; /* full relation (el caller verifica que sea válida) */
}

/*--------------------------------------------------------------------
 * try_combine_partial — Buscar match para una parcial y combinar.
 *
 * Para 1LP (lp2==0): buscar otra parcial con el mismo LP.
 * Para 2LP (lp2!=0): buscar otra parcial que comparta AL MENOS un LP.
 *   - Si comparte ambos: combinación directa → full relation.
 *   - Si comparte uno: la combinación tiene un LP residual.
 *     Podemos intentar combinar ESA con otra parcial (cadena).
 *     En esta versión, solo hacemos combinaciones simples (match directo).
 *
 * @return 0 = almacenada como parcial, 2 = combined → full relation
 *--------------------------------------------------------------------*/
int try_combine_partial(qs_struct *qs_data, mpz_t Qxi, mpz_t Xi,
                        int *exp_vec, int sign,
                        unsigned long lp1, unsigned long lp2)
{
    /* Buscar match en la tabla de parciales */
    long match_idx = -1;
    unsigned long shared_lp = 0;

    for (unsigned long k = 0; k < qs_data->partials.n; k++) {
        partial_entry *pe = &qs_data->partials.entries[k];

        if (lp2 == 0) {
            /* 1LP: buscar otra 1LP con el mismo primo grande */
            if (pe->large_prime2 == 0 && pe->large_prime == lp1) {
                match_idx = (long)k;
                shared_lp = lp1;
                break;
            }
        } else {
            /* 2LP: buscar otra parcial que comparta ambos LPs,
             * o al menos uno de ellos.
             * Prioridad 1: otra 2LP con los mismos dos LPs */
            if (pe->large_prime2 != 0) {
                if ((pe->large_prime == lp1 && pe->large_prime2 == lp2) ||
                    (pe->large_prime == lp2 && pe->large_prime2 == lp1)) {
                    /* Ambos LPs coinciden → full relation */
                    match_idx = (long)k;
                    shared_lp = lp1; /* ambos se cancelan */
                    break;
                }
            }
            /* Prioridad 2: otra 1LP con lp1 o lp2 → NO usamos esto
             * en la versión simple porque el resultado tendría un LP
             * residual y necesitaría otra ronda de combinación.
             * Lo dejamos para futuras mejoras. */
        }
    }

    if (match_idx >= 0) {
        partial_entry *match = &qs_data->partials.entries[match_idx];
        unsigned long remaining_lps[2];
        int n_remaining = 0;
        int rc = combine_two_partials(qs_data, exp_vec, sign, Qxi, Xi,
                                       match, shared_lp,
                                       remaining_lps, &n_remaining);

        /* Eliminar la parcial usada (swap con última) */
        free(match->exponents);
        mpz_clears(match->lhs, match->Qx, match->roota, match->a_value, NULL);
        unsigned long last = qs_data->partials.n - 1;
        if ((unsigned long)match_idx != last)
            qs_data->partials.entries[match_idx] = qs_data->partials.entries[last];
        qs_data->partials.n--;

        (void)rc;
        return 2; /* combined partial = full relation */
    } else {
        /* No hay match: guardar esta parcial */
        store_partial(qs_data, Qxi, Xi, exp_vec, sign, lp1, lp2);
        return 0;
    }
}

/*--------------------------------------------------------------------
 * classify_smooth_or_partial — decide qué hacer con un candidato cuyo
 * Qxi ya fue reducido a `residual` por la base.
 *
 * Toma ownership de exp_vec (lo libera para full/combined/discard;
 * lo transfiere a partial_entry si se almacena como 1LP/2LP no emparejada).
 * Toma ownership de residual (mpz_clear interno).
 *
 * @return 1 = full, 2 = combined partial = full, 0 = saved partial o discard.
 *--------------------------------------------------------------------*/
static int classify_smooth_or_partial(qs_struct *qs_data, mpz_t Qxi, mpz_t Xi,
                                       int *exp_vec, int sign, mpz_t residual)
{
    if (mpz_cmp_ui(residual, 1) == 0) {
        double _ti0 = omp_get_wtime();
        if (sign) insertarNumero(&qs_data->mat, qs_data->n_BSuaves, 0, 1);
        for (long i = 0; i < qs_data->base.length; i++)
            insertarNumero(&qs_data->mat, qs_data->n_BSuaves, i + 1, exp_vec[i] % 2);
        add_a_factors_to_matrix(qs_data);
        #pragma omp atomic update
        g_t_insert_matrix += omp_get_wtime() - _ti0;
        #pragma omp atomic update
        g_n_full++;
        free(exp_vec);
        mpz_clear(residual);
        return 1;
    }

    unsigned long lp1 = 0;
    if (mpz_fits_ulong_p(residual)) lp1 = mpz_get_ui(residual);

    int is_prime_1lp = 0;
    if (lp1 > 1 && lp1 < qs_data->large_prime_bound) {
        double _tp0 = omp_get_wtime();
        is_prime_1lp = (mpz_probab_prime_p(residual, 15) > 0);
        #pragma omp atomic update
        g_t_prime_test += omp_get_wtime() - _tp0;
    }
    if (is_prime_1lp) {
        #pragma omp atomic update
        g_n_1lp_try++;
        double _tc0 = omp_get_wtime();
        int rc = try_combine_partial(qs_data, Qxi, Xi, exp_vec, sign, lp1, 0);
        #pragma omp atomic update
        g_t_combine += omp_get_wtime() - _tc0;
        mpz_clear(residual);
        if (rc == 2) free(exp_vec);
        return rc;
    }

    if (qs_data->large_prime_bound2 > 0 && mpz_cmp_ui(residual, 1) > 0) {
        unsigned long long res_ull = 0;
        if (mpz_sizeinbase(residual, 2) <= 64) {
            if (mpz_fits_ulong_p(residual)) {
                res_ull = mpz_get_ui(residual);
            } else {
                mpz_t hi, lo;
                mpz_inits(hi, lo, NULL);
                mpz_tdiv_q_2exp(hi, residual, 32);
                mpz_tdiv_r_2exp(lo, residual, 32);
                res_ull = ((unsigned long long)mpz_get_ui(hi) << 32) |
                          (unsigned long long)mpz_get_ui(lo);
                mpz_clears(hi, lo, NULL);
            }
        }

        if (res_ull > qs_data->max_fb2 &&
            res_ull <= qs_data->large_prime_bound2 &&
            mpz_probab_prime_p(residual, 1) == 0) {
            #pragma omp atomic update
            g_n_2lp_try++;
            #pragma omp atomic update
            g_n_pollard_run++;
            unsigned long f1 = 0, f2 = 0;
            double _tp0 = omp_get_wtime();
            int found = factor_cofactor_pollard(residual, qs_data->large_prime_bound, &f1, &f2);
            #pragma omp atomic update
            g_t_pollard += omp_get_wtime() - _tp0;
            if (found && f1 > 1 && f2 > 1 &&
                f1 < qs_data->large_prime_bound &&
                f2 < qs_data->large_prime_bound) {
                if (f1 > f2) { unsigned long tmp = f1; f1 = f2; f2 = tmp; }
                double _tc0 = omp_get_wtime();
                int rc = try_combine_partial(qs_data, Qxi, Xi, exp_vec, sign, f1, f2);
                #pragma omp atomic update
                g_t_combine += omp_get_wtime() - _tc0;
                qs_data->n_dlp_stored++;
                if (rc == 2) qs_data->n_dlp_combined++;
                mpz_clear(residual);
                if (rc == 2) free(exp_vec);
                return rc;
            }
        }
    }

    mpz_clear(residual);
    free(exp_vec);
    return 0;
}

/*--------------------------------------------------------------------
 * trialDivisionRecip — Trial division usando recíprocos precomputados
 *
 * Test de divisibilidad: en vez de mpz_divisible_p (GMP genérico),
 * usamos el sieve_offset y las raíces de criba para determinar si
 * p divide Q(x). Si offset mod p == root1 o root2, entonces p | Q(x).
 *
 * El cálculo de offset mod p se hace con el recíproco precomputado:
 *   q = (uint32)(((uint64)(offset + rcorrect) * recip) >> 32)
 *   r = offset - q * p
 * Esto reemplaza una operación de remainder (~30 ciclos) por una
 * multiplicación + shift (~5 ciclos).
 *
 * Para la división real (cuando p sí divide), se usa mpz_tdiv_q_ui
 * que es la operación nativa de GMP para dividir por un unsigned long
 * (mucho más rápida que mpz_divexact con otro mpz_t).
 *
 * @param Qxi      Valor Q(x) a factorizar
 * @param qs_data  Estructura con base de primos y tabla de parciales
 * @param Xi       Valor x del candidato (para calcular lhs = a*x+b)
 * @param sieve_offset  Posición en el array de criba = x + xmax
 * @return 1 = full relation, 2 = combined partial, 0 = no relation
 *--------------------------------------------------------------------*/
int trialDivisionRecip(mpz_t Qxi, qs_struct *qs_data, mpz_t Xi,
                       unsigned long sieve_offset)
{
    double _ta, _tb, _td_acc = 0;
    double _t_func0 = omp_get_wtime();

    #pragma omp atomic update
    g_n_calls++;

    _ta = omp_get_wtime();
    int *exp_vec = (int *)calloc(qs_data->base.length, sizeof(int));
    mpz_t res;
    mpz_init(res);
    mpz_abs(res, Qxi);
    _tb = omp_get_wtime();
    #pragma omp atomic update
    g_t_calloc += _tb - _ta;

    int sign = (mpz_sgn(Qxi) < 0) ? 1 : 0;

    double _t_recip0 = omp_get_wtime();

    /* Tratar p=2 por separado: simplemente eliminar factores de 2 */
    if (qs_data->base.length > 0 && qs_data->base.primes[0].p == 2) {
        unsigned long twos = 0;
        double _td0 = omp_get_wtime();
        while (mpz_even_p(res)) {
            mpz_tdiv_q_2exp(res, res, 1);
            twos++;
        }
        _td_acc += omp_get_wtime() - _td0;
        exp_vec[0] = (int)twos;
    }

    /* Para cada primo de la base (excepto p=2): test con recíproco */
    long start_i = (qs_data->base.length > 0 && qs_data->base.primes[0].p == 2) ? 1 : 0;

    for (long i = start_i; i < qs_data->base.length; i++) {
        sieve_prime *sp = &qs_data->base.sp[i];
        uint32_t p     = sp->p;
        uint32_t root1 = sp->root1;
        uint32_t root2 = sp->root2;
        uint32_t recip = sp->recip;
        uint32_t rcorr = sp->rcorrect;

        if (root1 == UINT32_MAX) {
            /* Raíz inválida (p | a): hacer mod directo con GMP */
            if (mpz_divisible_ui_p(res, p)) {
                unsigned long cnt = 0;
                double _td0 = omp_get_wtime();
                do {
                    mpz_tdiv_q_ui(res, res, p);
                    cnt++;
                } while (mpz_divisible_ui_p(res, p));
                _td_acc += omp_get_wtime() - _td0;
                exp_vec[i] = (int)cnt;
            }
            continue;
        }

        /* Test de divisibilidad con recíproco:
         * q = (uint32)(((uint64)(sieve_offset + rcorr) * recip) >> 32)
         * remainder = sieve_offset - q * p
         * Si remainder == root1 o root2 → p divide Q(x) */
        uint32_t q = (uint32_t)(((uint64_t)(sieve_offset + rcorr) *
                                  (uint64_t)recip) >> 32);
        uint32_t remainder = (uint32_t)sieve_offset - q * p;

        if (remainder == root1 || remainder == root2) {
            /* p divide Q(x): hacer las divisiones sucesivas */
            unsigned long cnt = 0;
            double _td0 = omp_get_wtime();
            do {
                mpz_tdiv_q_ui(res, res, p);
                cnt++;
            } while (mpz_divisible_ui_p(res, p));
            _td_acc += omp_get_wtime() - _td0;
            exp_vec[i] = (int)cnt;
        }
        /* else: p no divide Q(x), skip (~5 ciclos) */
    }

    double _t_recip_total = omp_get_wtime() - _t_recip0;
    /* recip_loop excluye el tiempo de mpz_tdiv_q_ui, ya contabilizado en g_t_gmp_div */
    #pragma omp atomic update
    g_t_recip_loop += _t_recip_total - _td_acc;
    #pragma omp atomic update
    g_t_gmp_div    += _td_acc;

    int result = classify_smooth_or_partial(qs_data, Qxi, Xi, exp_vec, sign, res);
    #pragma omp atomic update
    g_t_trialDiv_total += omp_get_wtime() - _t_func0;
    return result;
}




/**
 * @brief Factoriza el array Qxi con divisiones usando recíprocos,
 * verifica si cada posicion es un numero bsuave y lo agrega al archivo polinomio.txt
 * @param qs_data estructura que contiene el array Qxi
 * @param endPos cantidad de candidatos
 * @param posXi índice inicial
 * @param xmax mitad del intervalo de criba (para calcular sieve_offset)
 * @return retorna 1 si aun faltan numeros B_suaves por verificar y 0 en caso de haberlos encontrado todos
 */
int factoringTrial(qs_struct * qs_data, unsigned long endPos, unsigned long posXi, unsigned long xmax){
	FILE *fp = g_polinomio_fp;
	for (unsigned long i = 0; i < endPos; i++)
	{
		long x_val = mpz_get_si(qs_data->intervalo.Xi[posXi]);
		unsigned long tf_offset = (unsigned long)((long)xmax + x_val);
		int result = trialDivisionRecip(qs_data->intervalo.Qxi[i], qs_data,
		                                qs_data->intervalo.Xi[posXi], tf_offset);
		if(result == 1){
			/* Full relation encontrada — vector ya insertado por trialDivision */
			qs_data->n_BSuaves++;
			/* Escribir lhs;a*Q(x);roota */
			double _tw0 = omp_get_wtime();
			mpz_t lhs, Qfile;
			mpz_inits(lhs, Qfile, NULL);
			mpz_mul(lhs, qs_data->poly.a, qs_data->intervalo.Xi[posXi]);
			mpz_add(lhs, lhs, qs_data->poly.b);
			/* Qfile = a * Q(x) = lhs² - N */
			mpz_mul(Qfile, lhs, lhs);
			mpz_sub(Qfile, Qfile, qs_data->n);
			mpz_out_str(fp, 10, lhs);
			fprintf(fp, ";");
			mpz_out_str(fp, 10, Qfile);
			fprintf(fp, ";");
			mpz_out_str(fp, 10, qs_data->roota);
			fprintf(fp, "\n");
			mpz_clears(lhs, Qfile, NULL);
			add_writeFullRel_time(omp_get_wtime() - _tw0);
			if(qs_data->n_BSuaves==qs_data->mat.n_rows){
				return 0;
			}
		} else if(result == 2){
			/* Combined partial — vector y polinomio.txt ya escritos por trialDivision */
			qs_data->n_BSuaves++;
			if(qs_data->n_BSuaves==qs_data->mat.n_rows){
				return 0;
			}
		}
		/* result == 0: no relation o parcial guardada */
		posXi++;
	}
	return 1;
}


