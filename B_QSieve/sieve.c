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

/* Tamaño del bloque de criba: 32KB cabe en L1 cache */
#define SIEVE_BLOCK_SIZE 32768

extern int CORES;

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

    for (long i = 0; i < qs_data->base.length; i++) {
        prime *fb = &qs_data->base.primes[i];
        fb->p = (uint32_t)mpz_get_ui(fb->value);
        /* logp = round(log2(p)) */
        if (fb->p >= 2)
            fb->logp = (uint8_t)(log2((double)fb->p) + 0.5);
        else
            fb->logp = 1;
        
        /* Precomputar sqrt(N) mod p */
        if (fb->p == 2) {
            fb->sqrt_n_mod_p = 1; /* N es impar, sqrt(N) mod 2 = 1 */
        } else {
            shanksTonelli(qs_data->n, fb->value, r1, r2);
            fb->sqrt_n_mod_p = (uint32_t)mpz_get_ui(r1);
        }
        fb->root1 = 0;
        fb->root2 = 0;
    }
    mpz_clears(r1, r2, NULL);
}

/*--------------------------------------------------------------------
 * sieve_compute_roots — Calcula raíces de criba para el polinomio actual.
 *
 * Para Q(x) = a*x² + 2*b*x + c, las raíces de Q(x) ≡ 0 (mod p) son:
 *   x ≡ a⁻¹ * (±sqrt(N) - b) (mod p)
 *
 * Las raíces se almacenan como offsets en [0, sieve_interval)
 * donde sieve_interval = 2 * xmax.
 *--------------------------------------------------------------------*/
void sieve_compute_roots(qs_struct *qs_data, unsigned long sieve_interval) {
    for (long i = 0; i < qs_data->base.length; i++) {
        prime *fb = &qs_data->base.primes[i];
        uint32_t p = fb->p;
        if (p < 2) {
            fb->root1 = fb->root2 = UINT32_MAX; /* inválido */
            continue;
        }

        /* a mod p */
        uint32_t a_mod_p = (uint32_t)mpz_fdiv_ui(qs_data->poly.a, p);
        uint32_t b_mod_p = (uint32_t)mpz_fdiv_ui(qs_data->poly.b, p);
        
        if (a_mod_p == 0) {
            /* p divide a → Q(x) es lineal mod p: 2*b*x + c ≡ 0 (mod p) */
            uint32_t twob = (2 * b_mod_p) % p;
            if (twob == 0) {
                fb->root1 = fb->root2 = UINT32_MAX;
                continue;
            }
            /* Inverso modular con extended gcd (uint32 nativo) */
            int64_t g, x0, y0;
            {
                int64_t aa = twob, bb = p, xx = 1, yy = 0, xx1 = 0, yy1 = 1, qq, tt;
                while (bb) { qq = aa/bb; tt = bb; bb = aa - qq*bb; aa = tt;
                    tt = xx1; xx1 = xx - qq*xx1; xx = tt;
                    tt = yy1; yy1 = yy - qq*yy1; yy = tt; }
                g = aa; x0 = xx; y0 = yy;
                (void)y0; (void)g;
            }
            uint32_t c_mod_p = (uint32_t)mpz_fdiv_ui(qs_data->poly.c, p);
            int64_t s = (-(int64_t)c_mod_p * x0) % (int64_t)p;
            if (s < 0) s += p;
            /* Ajustar al offset del intervalo: el intervalo va de -xmax a +xmax
             * posición en array = x + xmax, donde x es la raíz de criba.
             * La raíz s está en [0,p), necesitamos el primer offset ≥ 0 en el array.
             */
            unsigned long xmax = sieve_interval / 2;
            /* El primer x que satisface x ≡ s (mod p) y x >= -xmax es:
             * offset = (s + xmax) mod p */
            fb->root1 = (uint32_t)((s + xmax) % p);
            fb->root2 = UINT32_MAX; /* una sola raíz */
            continue;
        }
        
        /* a_inv mod p */
        uint32_t a_inv;
        {
            int64_t aa = a_mod_p, bb = p, xx = 1, xx1 = 0, qq, tt;
            while (bb) { qq = aa/bb; tt = bb; bb = aa - qq*bb; aa = tt;
                tt = xx1; xx1 = xx - qq*xx1; xx = tt; }
            a_inv = (uint32_t)(((xx % (int64_t)p) + p) % p);
        }
        
        uint32_t sqrt_n = fb->sqrt_n_mod_p;
        unsigned long xmax = sieve_interval / 2;
        
        /* root1 = a_inv * (sqrt_n - b) mod p */
        int64_t r1 = ((int64_t)sqrt_n - (int64_t)b_mod_p) % (int64_t)p;
        if (r1 < 0) r1 += p;
        r1 = ((int64_t)a_inv * r1) % p;
        
        /* root2 = a_inv * (-sqrt_n - b) mod p = a_inv * (p - sqrt_n - b) mod p */
        int64_t r2 = ((int64_t)(p - sqrt_n) - (int64_t)b_mod_p) % (int64_t)p;
        if (r2 < 0) r2 += p;
        r2 = ((int64_t)a_inv * r2) % p;
        
        /* Ajustar al offset en el array de criba [0, sieve_interval) */
        fb->root1 = (uint32_t)(((int64_t)r1 + (int64_t)xmax) % p);
        fb->root2 = (uint32_t)(((int64_t)r2 + (int64_t)xmax) % p);
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
    unsigned long sieve_interval = 2 * xmax; /* total de posiciones */
    unsigned long num_blocks = (sieve_interval + SIEVE_BLOCK_SIZE - 1) / SIEVE_BLOCK_SIZE;
    
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
        /* cutoff_config = 1.5 * error_bits, como en msieve para fb < 800 */
        unsigned int error_bits = 0;
        if (qs_data->large_prime_bound > 1) {
            error_bits = (unsigned int)(log2((double)qs_data->large_prime_bound) + 0.5);
        }
        unsigned int cutoff_config = (unsigned int)(1.5 * error_bits);
        
        /* bits(|c|) donde c = (b² - N) / a */
        mpz_t c_val, tmp;
        mpz_inits(c_val, tmp, NULL);
        mpz_mul(c_val, qs_data->poly.b, qs_data->poly.b);
        mpz_sub(c_val, c_val, qs_data->n);
        mpz_tdiv_q(c_val, c_val, qs_data->poly.a);
        mpz_abs(c_val, c_val);
        unsigned int c_bits = (unsigned int)mpz_sizeinbase(c_val, 2);
        mpz_clears(c_val, tmp, NULL);
        
        if (c_bits >= cutoff_config)
            cutoff = c_bits - cutoff_config;
        else
            cutoff = 0;
        
        /* Limitar a 250 (max uint8 útil) */
        if (cutoff > 250) cutoff = 250;
        if (cutoff < 2) cutoff = 2;
    }
    
    /* Calcular raíces de criba para este polinomio */
    sieve_compute_roots(qs_data, sieve_interval);
    
    /* Array temporal de criba: solo un bloque de 32KB en stack o malloc alineado */
    uint8_t *sieve_block = (uint8_t *)malloc(SIEVE_BLOCK_SIZE);
    if (!sieve_block) {
        fprintf(stderr, "Error al asignar sieve_block\n");
        exit(EXIT_FAILURE);
    }
    
    /* Buffer de candidatos */
    unsigned long capacity = sieve_interval / 20 + 256;
    long *indices = (long *)malloc(capacity * sizeof(long));
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
        
        /* Inicializar bloque con cutoff - 1 (como msieve) */
        memset(sieve_block, (uint8_t)(cutoff - 1), block_len);
        
        /* Cribar: restar logprime en cada posición de criba dentro del bloque */
        for (long i = 0; i < qs_data->base.length; i++) {
            prime *fb = &qs_data->base.primes[i];
            uint32_t p = fb->p;
            if (p < 2) continue;
            uint8_t logp = fb->logp;
            
            /* root1: primera posición dentro de este bloque */
            if (fb->root1 != UINT32_MAX) {
                uint32_t r1 = fb->root1;
                /* Avanzar r1 al rango [block_start, block_start + p) */
                if (r1 < block_start) {
                    unsigned long skip = (block_start - r1 + p - 1) / p;
                    r1 += (uint32_t)(skip * p);
                }
                for (uint32_t pos = r1; pos < block_end; pos += p) {
                    sieve_block[pos - block_start] -= logp;
                }
            }
            
            /* root2 */
            if (fb->root2 != UINT32_MAX && fb->root2 != fb->root1) {
                uint32_t r2 = fb->root2;
                if (r2 < block_start) {
                    unsigned long skip = (block_start - r2 + p - 1) / p;
                    r2 += (uint32_t)(skip * p);
                }
                for (uint32_t pos = r2; pos < block_end; pos += p) {
                    sieve_block[pos - block_start] -= logp;
                }
            }
        }
        
        /* Escanear el bloque: buscar posiciones con bit 7 set
         * (underflow = el valor original era >= cutoff y se restó bastante) */
        uint64_t *packed = (uint64_t *)sieve_block;
        unsigned long packed_len = block_len / 8;
        
        for (unsigned long qi = 0; qi < packed_len; qi++) {
            if ((packed[qi] & PACKED_MASK) == 0)
                continue;
            /* Hay al menos un candidato en estos 8 bytes */
            for (unsigned int jj = 0; jj < 8; jj++) {
                if (sieve_block[qi * 8 + jj] & 0x80) {
                    unsigned long pos = block_start + qi * 8 + jj;
                    if (pos < sieve_interval) {
                        long x = (long)pos - (long)xmax;
                        if (count >= capacity) {
                            capacity *= 2;
                            indices = (long *)realloc(indices, capacity * sizeof(long));
                        }
                        indices[count++] = x;
                    }
                }
            }
        }
        /* Bytes restantes (si block_len no es múltiplo de 8) */
        for (unsigned long qi = packed_len * 8; qi < block_len; qi++) {
            if (sieve_block[qi] & 0x80) {
                unsigned long pos = block_start + qi;
                if (pos < sieve_interval) {
                    long x = (long)pos - (long)xmax;
                    if (count >= capacity) {
                        capacity *= 2;
                        indices = (long *)realloc(indices, capacity * sizeof(long));
                    }
                    indices[count++] = x;
                }
            }
        }
    }
    
    free(sieve_block);
    *out_indices = indices;
    *out_count = count;
}
