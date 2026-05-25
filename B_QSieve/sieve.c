/*--------------------------------------------------------------------
 * sieve.c — Criba logarítmica uint8 cache-friendly (32KB blocks)
 *
 * Técnica msieve: array uint8 inicializado a cutoff, se resta logprime
 * en cada posición de criba. Las posiciones cuyo bit 7 queda encendido
 * (underflow = valor suave) son candidatas para trial division.
 *
 * Las raíces de criba se precomputan como uint32 nativos y se actualizan
 * incrementalmente entre polinomios derivados (Gray code).
 *--------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <mpfr.h>
#include "structsqs.h"
#include <string.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#ifdef __AVX2__
#include <immintrin.h>
#endif

/* Tamaño del bloque de criba: 32KB cabe en L1 cache */
#define SIEVE_BLOCK_SIZE 32768

/* Primos p <= SMALL_PRIME_MAX no se criban: golpean una fracción enorme de
 * posiciones (p=2 toca la mitad) pero aportan poco log. Se omiten del bucle
 * interno y se compensa bajando el umbral. Trial division sigue dividiendo
 * por ellos, así que la corrección de relaciones no se ve afectada. */
#define SMALL_PRIME_MAX 31
/* Multiplicador de la corrección de umbral (ver g_small_correction). Saltar
 * primos medianos vuelve la criba imprecisa (explosión de candidatos), por
 * eso el cutoff se queda chico — solo los p<=31, como msieve. */
#define SMALL_PRIME_CORR 1.3

/* Entry compacta para el inner loop de criba: p, logp y dos posiciones con
 * carry-forward entre bloques (evita la división por bloque). 16 bytes. */
typedef struct {
    uint32_t p;
    uint32_t nloc1;   /* posición de root1 dentro del bloque actual (carry) */
    uint32_t nloc2;   /* posición de root2, o UINT32_MAX = raíz única/inválida */
    uint8_t  logp;
    uint8_t  _pad[3];
} packed_sieve_t;

extern int CORES;

/* Calculado una vez: nº de primos iniciales con p <= SMALL_PRIME_MAX y la
 * corrección de log esperada que aportarían a una posición suave. */
static long g_sieve_skip = -1;
static unsigned int g_small_correction = 0;

/* Caches persistentes entre llamadas a sieve_mpqs (una factorización = tamaño fijo). */
static packed_sieve_t *g_psieve       = NULL;
static long            g_psieve_cap   = 0;
static uint8_t        *g_sieve_block  = NULL;
static long           *g_cand_buf     = NULL; /* pool de candidatos — sin malloc por llamada */
static unsigned long   g_cand_cap     = 0;

/* Acumuladores para diagnóstico del bottleneck del sieve */
static double g_t_sieve_cutoff = 0;
static double g_t_sieve_roots  = 0;
static double g_t_sieve_loop   = 0;
static double g_t_sieve_scan   = 0;
static unsigned long g_n_sieve_calls = 0;

void print_sieve_stats(double total_wall) {
    if (total_wall <= 0) total_wall = 1e-9;
    double sum = g_t_sieve_cutoff + g_t_sieve_roots + g_t_sieve_loop + g_t_sieve_scan;
    double other = total_wall - sum;
    if (other < 0) other = 0;
    printf("\n  [sieve_mpqs breakdown — %lu calls]\n", g_n_sieve_calls);
    printf("    cutoff calc          : %.3fs (%5.1f%%)\n", g_t_sieve_cutoff, 100.0*g_t_sieve_cutoff/total_wall);
    printf("    sieve_compute_roots  : %.3fs (%5.1f%%)\n", g_t_sieve_roots,  100.0*g_t_sieve_roots/total_wall);
    printf("    memset+sieve loop    : %.3fs (%5.1f%%)\n", g_t_sieve_loop,   100.0*g_t_sieve_loop/total_wall);
    printf("    scan candidates      : %.3fs (%5.1f%%)\n", g_t_sieve_scan,   100.0*g_t_sieve_scan/total_wall);
    printf("    otros (malloc/free)  : %.3fs (%5.1f%%)\n", other,            100.0*other/total_wall);
    fflush(stdout);
}

/*--------------------------------------------------------------------
 * shanksTonelli — raíz cuadrada modular (sin cambios)
 *--------------------------------------------------------------------*/
int shanksTonelli(mpz_t n, mpz_t p, mpz_t r1, mpz_t r2) {
	
	mpz_t resMod, p1, div;
	mpz_inits(resMod,p1,div,NULL);
	
	mpz_t q, ss;
	mpz_inits(q,ss,NULL);
	mpz_sub_ui(q,p,1);
	
	while (mpz_divisible_ui_p(q,2) != 0)
	{
		mpz_add_ui(ss,ss,1);
		mpz_divexact_ui(q,q,2);
	}

	if (mpz_cmp_ui(ss,1) == 0)
	{
		mpz_add_ui(p1,p,1);
		mpz_divexact_ui(div,p1,4);
		mpz_powm(r1,n,div,p);
		mpz_sub(r2,p,r1);
		
	}else{
		mpz_sub_ui(p1,p,1);
		
		mpz_t z;
		mpz_init(z);
		mpz_set_ui(z,2);
		
		mpz_divexact_ui(div,p1,2);
		
		mpz_powm(resMod,z,div,p);
		
		while (mpz_cmp(resMod,p1)!=0)
		{
			mpz_add_ui(z,z,1);
			mpz_powm(resMod,z,div,p);
		}
		
		mpz_t c,r,t,m;
		mpz_inits(c,r,t,m,NULL);
		
		mpz_powm(c,z,q,p);
		mpz_powm(t,n,q,p);
		
		mpz_add_ui(q,q,1);
		mpz_divexact_ui(div,q,2);
		
		mpz_powm(r,n,div,p);
		
		mpz_set(m,ss);
		
		while (1)
		{
			if (mpz_cmp_ui(t,1) == 0)
			{
				mpz_set(r1,r);
				mpz_t pr;
				mpz_init(pr);
				mpz_sub(pr,p,r);
				mpz_set(r2,pr);
				mpz_clear(pr);
				break;
			}
			
			mpz_t i,zz,m1;
			mpz_inits(i,zz,m1,NULL);
			
			mpz_set_ui(i,0);
			mpz_set(zz,t);
			
			mpz_sub_ui(m1,m,1);
			
			while (mpz_cmp_ui(zz,1) != 0 && mpz_cmp(i,m1) < 0)
			{
				mpz_powm_ui(zz,zz,2,p);
				mpz_add_ui(i,i,1);
			}
			
			mpz_t b,e;
			mpz_inits(b,e,NULL);
			
			mpz_set(b,c);
			mpz_set(e,m);
			mpz_sub(e,e,i);
			mpz_sub_ui(e,e,1);
			
			while (mpz_cmp_ui(e,0) > 0)
			{
				mpz_powm_ui(b,b,2,p);
				mpz_sub_ui(e,e,1);
			}
			
			mpz_mul(r,r,b);
			mpz_powm_ui(r,r,1,p);
			
			mpz_powm_ui(c,b,2,p);
			
			mpz_mul(t,t,c);
			mpz_powm_ui(t,t,1,p);

			mpz_set(m,i);
			mpz_clears(i,zz,m1,b,e,NULL);
		}
		mpz_clears(c,r,t,m,z,NULL);
	}
	mpz_clears(resMod,p1,div,q,ss,NULL);
	return 1;
}

/*--------------------------------------------------------------------
 * sieve_precompute_roots — Precomputa sqrt(N) mod p y campos nativos
 * para cada primo de la base. Se llama UNA VEZ tras generar la base.
 *--------------------------------------------------------------------*/
void sieve_precompute_roots(qs_struct *qs_data) {
    mpz_t r1, r2;
    mpz_inits(r1, r2, NULL);

    /* Asignar el array compacto paralelo a la base */
    qs_data->base.sp = (sieve_prime *)malloc(
        (size_t)qs_data->base.length * sizeof(sieve_prime));

    for (long i = 0; i < qs_data->base.length; i++) {
        prime *fb = &qs_data->base.primes[i];
        sieve_prime *sp = &qs_data->base.sp[i];

        uint32_t p = (uint32_t)mpz_get_ui(fb->value);
        fb->p = p;
        sp->p = p;
        sp->logp = (p >= 2) ? (uint8_t)(log2((double)p) + 0.5) : 1;

        /* Precomputar sqrt(N) mod p */
        if (p == 2) {
            fb->sqrt_n_mod_p = 1; /* N es impar, sqrt(N) mod 2 = 1 */
        } else {
            shanksTonelli(qs_data->n, fb->value, r1, r2);
            fb->sqrt_n_mod_p = (uint32_t)mpz_get_ui(r1);
        }
        sp->root1 = 0;
        sp->root2 = 0;

        /* Recíproco para trial division: recip = ⌊2^32/p⌋ (o +1 si el
         * truncamiento pierde precisión). Permite calcular (x % p) como
         *   q = (uint32)((uint64)x * recip >> 32);  r = x - q*p;
         * Ver Agner Fog, "Optimizing subroutines in assembly language". */
        if (p >= 2) {
            uint64_t r64 = ((uint64_t)1 << 32) / (uint64_t)p;
            double exact = 4294967296.0 / (double)p;  /* 2^32 / p */
            if (fabs(exact - (double)r64) < 0.5) {
                sp->rcorrect = 1;  /* recip exacto (o redondeado abajo) */
                sp->recip = (uint32_t)r64;
            } else {
                sp->rcorrect = 0;  /* necesita +1 para corregir */
                sp->recip = (uint32_t)(r64 + 1);
            }
        } else {
            sp->recip = 0;
            sp->rcorrect = 0;
        }
    }
    mpz_clears(r1, r2, NULL);
}

/* Deltas Gray code: g_root_delta[j*base.length + i] = a⁻¹·(2·B_j) mod p_i.
 * Permite actualizar las raíces de criba entre polinomios derivados con un
 * solo add+mod por primo en vez de recomputar a⁻¹ y las raíces desde cero. */
static uint32_t *g_root_delta = NULL;
static long g_root_delta_cap = 0;

/* Inverso modular de v mod p vía extended gcd (uint32 nativo). */
static inline uint32_t modinv_u32(uint32_t v, uint32_t p) {
    int64_t aa = v, bb = p, xx = 1, xx1 = 0, qq, tt;
    while (bb) { qq = aa/bb; tt = bb; bb = aa - qq*bb; aa = tt;
        tt = xx1; xx1 = xx - qq*xx1; xx = tt; }
    return (uint32_t)(((xx % (int64_t)p) + p) % p);
}

/* Raíz de criba para un primo q que divide 'a' (Q(x) ≡ 2bx+c mod q es lineal).
 * Devuelve el offset en [0,p) o UINT32_MAX si 2b ≡ 0 (mod q). */
static uint32_t linear_sieve_root(qs_struct *qs_data, uint32_t p, unsigned long xmax) {
    uint32_t b_mod_p = (uint32_t)mpz_fdiv_ui(qs_data->poly.b, p);
    uint32_t twob = (2 * b_mod_p) % p;
    if (twob == 0) return UINT32_MAX;
    uint32_t twob_inv = modinv_u32(twob, p);
    uint32_t c_mod_p = (uint32_t)mpz_fdiv_ui(qs_data->poly.c, p);
    int64_t s = (-(int64_t)c_mod_p * (int64_t)twob_inv) % (int64_t)p;
    if (s < 0) s += p;
    return (uint32_t)(((unsigned long)s + xmax) % p);
}

/*--------------------------------------------------------------------
 * sieve_init_roots_for_a — Cómputo completo de raíces + precompute de
 * deltas Gray code. Se llama una vez por cada nuevo 'a'.
 *
 * Para Q(x) = a*x² + 2*b*x + c, las raíces de Q(x) ≡ 0 (mod p) son
 *   x ≡ a⁻¹·(±sqrt(N) - b) (mod p),
 * almacenadas como offsets en [0,p) tras sumar xmax.
 *--------------------------------------------------------------------*/
void sieve_init_roots_for_a(qs_struct *qs_data, unsigned long sieve_interval) {
    long blen = qs_data->base.length;
    unsigned int s = qs_data->siqs_state.num_factors;
    unsigned long xmax = sieve_interval / 2;

    if (g_root_delta_cap < blen) {
        free(g_root_delta);
        g_root_delta = (uint32_t *)malloc((size_t)blen * MAX_SIQS_FACTORS * sizeof(uint32_t));
        g_root_delta_cap = blen;
    }

    for (long i = 0; i < blen; i++) {
        prime *fb = &qs_data->base.primes[i];
        sieve_prime *sp = &qs_data->base.sp[i];
        uint32_t p = sp->p;
        if (p < 2) {
            sp->root1 = sp->root2 = UINT32_MAX;
            for (unsigned int j = 0; j < s; j++)
                g_root_delta[(size_t)j * blen + i] = 0;
            continue;
        }

        uint32_t a_mod_p = (uint32_t)mpz_fdiv_ui(qs_data->poly.a, p);

        if (a_mod_p == 0) {
            /* p | a → caso lineal; sin deltas (se recomputa cada poly) */
            sp->root1 = linear_sieve_root(qs_data, p, xmax);
            sp->root2 = UINT32_MAX;
            for (unsigned int j = 0; j < s; j++)
                g_root_delta[(size_t)j * blen + i] = 0;
            continue;
        }

        uint32_t a_inv = modinv_u32(a_mod_p, p);
        uint32_t b_mod_p = (uint32_t)mpz_fdiv_ui(qs_data->poly.b, p);
        uint32_t sqrt_n = fb->sqrt_n_mod_p;

        int64_t r1 = ((int64_t)sqrt_n - (int64_t)b_mod_p) % (int64_t)p;
        if (r1 < 0) r1 += p;
        r1 = ((int64_t)a_inv * r1) % p;

        int64_t r2 = ((int64_t)(p - sqrt_n) - (int64_t)b_mod_p) % (int64_t)p;
        if (r2 < 0) r2 += p;
        r2 = ((int64_t)a_inv * r2) % p;

        sp->root1 = (uint32_t)((r1 + (int64_t)xmax) % p);
        sp->root2 = (uint32_t)((r2 + (int64_t)xmax) % p);

        /* delta[j] = a⁻¹·(2·B_j) mod p — B_j ya viene doblado en Bvals */
        for (unsigned int j = 0; j < s; j++) {
            uint32_t bj = (uint32_t)mpz_fdiv_ui(qs_data->siqs_state.Bvals[j], p);
            g_root_delta[(size_t)j * blen + i] =
                (uint32_t)(((uint64_t)a_inv * bj) % p);
        }
    }
}

/*--------------------------------------------------------------------
 * sieve_update_roots — Actualización incremental de raíces entre
 * polinomios derivados. b cambió por sign·(2·B_j), así que cada raíz
 * cambia por -sign·delta[j] (mod p): un add + resta condicional por primo.
 *--------------------------------------------------------------------*/
void sieve_update_roots(qs_struct *qs_data, unsigned long sieve_interval) {
    long blen = qs_data->base.length;
    unsigned int j = qs_data->siqs_state.sieve_flip_j;
    int sign = qs_data->siqs_state.sieve_flip_sign;
    const uint32_t *delta_j = &g_root_delta[(size_t)j * blen];

    for (long i = 0; i < blen; i++) {
        sieve_prime *sp = &qs_data->base.sp[i];
        uint32_t p = sp->p;
        if (p < 2) continue;
        uint32_t d = delta_j[i];
        uint32_t shift = (sign > 0) ? (p - d) : d; /* root += shift (mod p) */
        if (sp->root1 != UINT32_MAX) {
            uint32_t r = sp->root1 + shift;
            if (r >= p) r -= p;
            sp->root1 = r;
        }
        if (sp->root2 != UINT32_MAX) {
            uint32_t r = sp->root2 + shift;
            if (r >= p) r -= p;
            sp->root2 = r;
        }
    }

    /* Los s primos que dividen 'a' (caso lineal) sí se recomputan: su raíz
     * depende de b y c de forma no incremental. Son pocos (s ≤ 9). */
    unsigned long xmax = sieve_interval / 2;
    unsigned int s = qs_data->siqs_state.num_factors;
    for (unsigned int k = 0; k < s; k++) {
        long idx = (long)qs_data->siqs_state.factor_fb_idx[k];
        sieve_prime *sp = &qs_data->base.sp[idx];
        sp->root1 = linear_sieve_root(qs_data, sp->p, xmax);
        sp->root2 = UINT32_MAX;
    }
}

/*--------------------------------------------------------------------
 * sieve_mpqs — Criba logarítmica uint8 en bloques de 32KB
 *
 * Técnica msieve: inicializar el bloque con un valor de cutoff, luego
 * restar logprime en cada posición de criba. Las posiciones donde el
 * byte "underflows" (bit 7 set) son candidatas a ser suaves.
 *
 * El intervalo total es sieve_interval = 2*xmax posiciones.
 * Se divide en bloques de SIEVE_BLOCK_SIZE bytes.
 *--------------------------------------------------------------------*/
void sieve_mpqs(qs_struct *qs_data, unsigned long xmax,
                long **out_indices, unsigned long *out_count)
{
    g_n_sieve_calls++;
    double _t0 = omp_get_wtime();

    unsigned long sieve_interval = 2 * xmax; /* total de posiciones */
    unsigned long num_blocks = (sieve_interval + SIEVE_BLOCK_SIZE - 1) / SIEVE_BLOCK_SIZE;

    /* Primera llamada: contar primos chicos a saltar y su corrección de log.
     * Una posición es golpeada por p con prob ~n_roots/p y al serlo resta
     * logp; la corrección es la suma esperada de esos aportes. */
    if (g_sieve_skip < 0) {
        g_sieve_skip = 0;
        double corr = 0.0;
        for (long i = 0; i < qs_data->base.length; i++) {
            uint32_t p = qs_data->base.sp[i].p;
            if (p > SMALL_PRIME_MAX) break;
            g_sieve_skip++;
            int nroots = (p == 2) ? 1 : 2;
            corr += (double)nroots / (double)p * (double)qs_data->base.sp[i].logp;
        }
        /* Factor SMALL_PRIME_CORR: los números suaves son más divisibles por
         * primos chicos que el promedio, así que pierden más log que el valor
         * esperado. Se re-tunea según SMALL_PRIME_MAX. */
        g_small_correction = (unsigned int)(corr * SMALL_PRIME_CORR + 0.5);
    }

    /* Calcular cutoff al estilo msieve:
     * 
     * El sieve array se inicializa con (cutoff_fill - 1). Se resta log2(p)
     * por cada primo que divide la posición. Si la suma de logs supera
     * cutoff_fill, hay underflow → bit 7 set → candidato.
     *
     * cutoff_fill = bits(|c|) - cutoff_config
     * donde c = (b² - N)/a es el coeficiente constante del polinomio
     * y cutoff_config = 1.5 * log2(LP_bound) para bases pequeñas (<800 primos)
     *
     * Esto acepta candidatos donde el cofactor residual es < LP_bound^1.5
     * (generoso para capturar relaciones con 1 primo grande).
     */
    unsigned int cutoff;
    {
        /* cutoff_config precomputado en main (= 1.5 * log2(LP_bound)) */
        unsigned int cutoff_config = qs_data->sieve_cutoff_config;

        /* bits(|c|) donde c = (b² - N) / a */
        mpz_t c_val;
        mpz_init(c_val);
        mpz_mul(c_val, qs_data->poly.b, qs_data->poly.b);
        mpz_sub(c_val, c_val, qs_data->n);
        mpz_tdiv_q(c_val, c_val, qs_data->poly.a);
        mpz_abs(c_val, c_val);
        unsigned int c_bits = (unsigned int)mpz_sizeinbase(c_val, 2);
        mpz_clear(c_val);
        
        if (c_bits >= cutoff_config)
            cutoff = c_bits - cutoff_config;
        else
            cutoff = 0;

        /* Compensar los primos chicos no cribados: una posición suave acumula
         * g_small_correction menos log, así que se baja el umbral igual. */
        if (cutoff > g_small_correction)
            cutoff -= g_small_correction;

        /* Limitar a 250 (max uint8 útil) */
        if (cutoff > 250) cutoff = 250;
        if (cutoff < 2) cutoff = 2;
    }
    double _t1 = omp_get_wtime();
    g_t_sieve_cutoff += _t1 - _t0;

    /* Raíces de criba: cómputo completo si cambió 'a', incremental (Gray
     * code) si solo cambió 'b' respecto al polinomio derivado anterior. */
    if (qs_data->siqs_state.sieve_new_a)
        sieve_init_roots_for_a(qs_data, sieve_interval);
    else
        sieve_update_roots(qs_data, sieve_interval);
    double _t2 = omp_get_wtime();
    g_t_sieve_roots += _t2 - _t1;

    /* Array compacto con carry-forward de raíces: evita recalcular el offset
     * por división al inicio de cada bloque. Reutilizado entre llamadas;
     * solo se reasigna si el FB creció (nunca ocurre dentro de una factorización). */
    long n_sp = qs_data->base.length - g_sieve_skip;
    if (n_sp > g_psieve_cap) {
        free(g_psieve);
        g_psieve     = (packed_sieve_t *)malloc(n_sp * sizeof(packed_sieve_t));
        g_psieve_cap = n_sp;
        if (!g_psieve) { fprintf(stderr, "Error al asignar psieve\n"); exit(EXIT_FAILURE); }
    }
    packed_sieve_t *psieve = g_psieve;
    {
        sieve_prime *sp = qs_data->base.sp + g_sieve_skip;
        for (long i = 0; i < n_sp; i++, sp++) {
            uint32_t r1 = sp->root1;
            uint32_t r2 = sp->root2;
            /* Ordenar r1 ≤ r2; UINT32_MAX señaliza raíz inválida/única */
            if (r2 != UINT32_MAX && r2 < r1) { uint32_t t = r1; r1 = r2; r2 = t; }
            psieve[i].p     = sp->p;
            psieve[i].logp  = sp->logp;
            psieve[i].nloc1 = r1;
            psieve[i].nloc2 = r2;
        }
    }

    /* Bloque de criba: asignado una sola vez, reutilizado entre llamadas. */
    if (!g_sieve_block) {
        g_sieve_block = (uint8_t *)malloc(SIEVE_BLOCK_SIZE);
        if (!g_sieve_block) { fprintf(stderr, "Error al asignar sieve_block\n"); exit(EXIT_FAILURE); }
    }
    uint8_t *sieve_block = g_sieve_block;
    
    /* Buffer de candidatos — pool persistente, sin malloc por llamada */
    unsigned long capacity = sieve_interval / 20 + 256;
    if (capacity > g_cand_cap) {
        free(g_cand_buf);
        g_cand_buf = (long *)malloc(capacity * sizeof(long));
        g_cand_cap = capacity;
    }
    long *indices = g_cand_buf;
    unsigned long count = 0;
    
    /* Macro para detectar si algún byte de un uint64 tiene bit 7 set */
    #define PACKED_MASK 0x8080808080808080ULL
    
    /* Para cada bloque de criba */
    for (unsigned long blk = 0; blk < num_blocks; blk++) {
        unsigned long block_start = blk * SIEVE_BLOCK_SIZE;
        unsigned long block_end = block_start + SIEVE_BLOCK_SIZE;
        if (block_end > sieve_interval)
            block_end = sieve_interval;
        unsigned long block_len = block_end - block_start;

        double _tl0 = omp_get_wtime();
        /* Inicializar bloque completo (incluso bytes más allá del intervalo real
         * en el último bloque parcial: se criban pero nunca se escanean). */
        memset(sieve_block, (uint8_t)(cutoff - 1), SIEVE_BLOCK_SIZE);

        /* Cribar con carry-forward: psieve[i].nloc1/nloc2 ya apuntan al inicio
         * correcto dentro de este bloque — sin división.
         *   Two-root path: loop fusionado root1+root2 al estilo msieve.
         *   Single-root path: loop simple (primos que dividen 'a', ≤ s ≤ 9). */
        for (long i = 0; i < n_sp; i++) {
            uint32_t p    = psieve[i].p;
            uint8_t  logp = psieve[i].logp;
            uint32_t r1   = psieve[i].nloc1;
            uint32_t r2   = psieve[i].nloc2;

            if (r2 == UINT32_MAX) {
                /* raíz única: loop simple */
                while (r1 < SIEVE_BLOCK_SIZE) {
                    sieve_block[r1] -= logp;
                    r1 += p;
                }
                psieve[i].nloc1 = r1 - SIEVE_BLOCK_SIZE;
            } else if (p >= (SIEVE_BLOCK_SIZE / 2)) {
                /* primo grande: cada raíz golpea ≤2 veces por bloque */
                while (r1 < SIEVE_BLOCK_SIZE) { sieve_block[r1] -= logp; r1 += p; }
                while (r2 < SIEVE_BLOCK_SIZE) { sieve_block[r2] -= logp; r2 += p; }
                psieve[i].nloc1 = r1 - SIEVE_BLOCK_SIZE;
                psieve[i].nloc2 = r2 - SIEVE_BLOCK_SIZE;
            } else {
                /* primo medio: dos raíces, unrolled ×2.
                 * Re-sort r1 ≤ r2: el residual del bloque anterior puede romper
                 * el invariante. La condición r2+p < BLOCK garantiza r1+p < BLOCK. */
                if (r1 > r2) { uint32_t _t = r1; r1 = r2; r2 = _t; }
                uint32_t p2 = p + p;
                while (r2 + p < SIEVE_BLOCK_SIZE) {
                    sieve_block[r1]     -= logp;
                    sieve_block[r2]     -= logp;
                    sieve_block[r1 + p] -= logp;
                    sieve_block[r2 + p] -= logp;
                    r1 += p2;
                    r2 += p2;
                }
                /* residual independiente por raíz (≤1 iter cada una) */
                while (r1 < SIEVE_BLOCK_SIZE) { sieve_block[r1] -= logp; r1 += p; }
                while (r2 < SIEVE_BLOCK_SIZE) { sieve_block[r2] -= logp; r2 += p; }
                psieve[i].nloc1 = r1 - SIEVE_BLOCK_SIZE;
                psieve[i].nloc2 = r2 - SIEVE_BLOCK_SIZE;
            }
        }

        double _tl1 = omp_get_wtime();
        g_t_sieve_loop += _tl1 - _tl0;

        /* Escanear el bloque: buscar posiciones con bit 7 set.
         * AVX2: _mm256_movemask_epi8 extrae MSB de 32 bytes → máscara de 32 bits,
         *       __builtin_ctz localiza cada bit encendido en O(hits) en vez de O(32).
         * pos < sieve_interval es siempre cierto dentro de [0, block_len). */
#ifdef __AVX2__
        {
            unsigned long avx_steps = block_len / 32;
            for (unsigned long qi = 0; qi < avx_steps; qi++) {
                uint32_t mask = (uint32_t)_mm256_movemask_epi8(
                    _mm256_loadu_si256((const __m256i *)(sieve_block + qi * 32)));
                while (mask) {
                    int bit = __builtin_ctz(mask);
                    indices[count++] = (long)(block_start + qi * 32 + (unsigned)bit) - (long)xmax;
                    mask &= mask - 1;
                }
            }
            for (unsigned long qi = avx_steps * 32; qi < block_len; qi++) {
                if (sieve_block[qi] & 0x80)
                    indices[count++] = (long)(block_start + qi) - (long)xmax;
            }
        }
#else
        {
            uint64_t *packed = (uint64_t *)sieve_block;
            unsigned long packed_len = block_len / 8;
            for (unsigned long qi = 0; qi < packed_len; qi++) {
                if ((packed[qi] & PACKED_MASK) == 0) continue;
                for (unsigned int jj = 0; jj < 8; jj++) {
                    if (sieve_block[qi * 8 + jj] & 0x80)
                        indices[count++] = (long)(block_start + qi * 8 + jj) - (long)xmax;
                }
            }
            for (unsigned long qi = packed_len * 8; qi < block_len; qi++) {
                if (sieve_block[qi] & 0x80)
                    indices[count++] = (long)(block_start + qi) - (long)xmax;
            }
        }
#endif
        g_t_sieve_scan += omp_get_wtime() - _tl1;
    }

    *out_indices = indices;
    *out_count = count;
}
