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

//Cantidad de procesadores logicos que se quieren usar definidos en B_QSieve
extern int CORES;

/** 
* @brief 
* @param n: 
* @param p: numero primo 
* @param r1: root 1 parametro de salida
* @param r2: root 2 parametro de salida
* @return 1 si se encontro solucion o 0 si no se encuentra una solucion
*/
int shanksTonelli(mpz_t n, mpz_t p, mpz_t r1, mpz_t r2) {
	
	mpz_t resMod, p1, div;//p1 sera prime-1, resMod y div seran variables temporales para resultados de operaciones
	mpz_inits(resMod,p1,div,NULL);
	
	/*Se puede omitir esta parte por que los primos que llegan 
	 * ya se les calculo el simbolo de legendre
	if((mpz_legendre(n,p)==1)){
		mpz_set_ui(r1,0);
		mpz_set_ui(r2,0);
		return 0;
	}*/
	
	mpz_t q, ss;
	mpz_inits(q,ss,NULL);
	mpz_sub_ui(q,p,1);//q=(prime-1)
	
	
	//mientras que el ultimo bit de q sea 0 (q par)
	while (mpz_divisible_ui_p(q,2) != 0)
	{
		mpz_add_ui(ss,ss,1);
		mpz_divexact_ui(q,q,2);
	}

	if (mpz_cmp_ui(ss,1) == 0)
	{
		mpz_add_ui(p1,p,1);//p1=prime+1;
		mpz_divexact_ui(div,p1,4);
		//div=(prime+1)/4
		mpz_powm(r1,n,div,p);//r1=n^((prime+1)/4) mod prime
		mpz_sub(r2,p,r1);//r2=prime-r1
		
	}else{
		mpz_sub_ui(p1,p,1);//p=prime-1;
		
		mpz_t z;
		mpz_init(z);
		mpz_set_ui(z,2);
		
		mpz_divexact_ui(div,p1,2);
		//(prime-1)/2
		
		mpz_powm(resMod,z,div,p);//resMod=z^(div) mod prime
		
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
		mpz_divexact_ui(div,q,2);//div=q/2
		
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
				mpz_set(r2,pr);//r2=prime-r
				mpz_clear(pr);
				break;
			}
			
			mpz_t i,zz,m1;
			mpz_inits(i,zz,m1,NULL);
			
			mpz_set_ui(i,0);
			mpz_set(zz,t);
			
			mpz_sub_ui(m1,m,1);
			
			//zz != 1 && i < (m-1)
			while (mpz_cmp_ui(zz,1) != 0 && mpz_cmp(i,m1) < 0)
			{
				mpz_powm_ui(zz,zz,2,p);//zz=zz*zz mod prime
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
				mpz_powm_ui(b,b,2,p);//b=b*b mod prime
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

/**
 * @brief Criba logarítmica optimizada para un polinomio MPQS.
 *
 * Para el polinomio Q(x) = a*x² + 2*b*x + c, calcula las dos raíces
 * de criba s₁, s₂ para cada primo p de la base (soluciones de Q(x) ≡ 0 mod p),
 * y acumula log(p) en un array de tamaño 2*xmax+1 en las posiciones
 * s₁, s₁+p, s₁+2p, … y s₂, s₂+p, s₂+2p, …
 *
 * Solo los índices cuya suma de logs supera el umbral T se devuelven
 * como candidatos para trial division.
 *
 * @param qs_data       Estructura con base de primos y polinomio MPQS
 * @param xmax          Mitad del intervalo [-xmax .. +xmax]
 * @param out_indices   Array de salida con las posiciones x candidatas (caller libera)
 * @param out_count     Número de candidatos
 */
void sieve_mpqs(qs_struct *qs_data, unsigned long xmax,
                long **out_indices, unsigned long *out_count)
{
    unsigned long total = 2 * xmax + 1;  // índices 0..total-1 representan x = -xmax..+xmax

    // Array de criba: S[i] acumula sum(log(p)) para x = i - xmax
    float *S = (float *)calloc(total, sizeof(float));
    if (!S) {
        fprintf(stderr, "Error al asignar memoria para array de criba\n");
        exit(EXIT_FAILURE);
    }

    // Calcular umbral T = ln(sqrt(|c|) * xmax) con margen
    // Heurística: T ≈ ln(sqrt(max|Q(x)|)) - ln(p_max)
    // max|Q(x)| ≈ a*xmax² + 2*|b|*xmax + |c|
    mpfr_t T_thr;
    mpfr_init2(T_thr, 128);
    {
        mpz_t maxQ;
        mpz_init(maxQ);
        // maxQ = a * xmax^2
        mpz_set_ui(maxQ, xmax);
        mpz_mul_ui(maxQ, maxQ, xmax);
        mpz_mul(maxQ, maxQ, qs_data->poly.a);
        // + 2*|b|*xmax
        mpz_t tmp;
        mpz_init(tmp);
        mpz_abs(tmp, qs_data->poly.b);
        mpz_mul_ui(tmp, tmp, 2 * xmax);
        mpz_add(maxQ, maxQ, tmp);
        // + |c|
        mpz_abs(tmp, qs_data->poly.c);
        mpz_add(maxQ, maxQ, tmp);

        mpfr_set_z(T_thr, maxQ, MPFR_RNDN);
        mpfr_sqrt(T_thr, T_thr, MPFR_RNDZ);
        mpfr_log(T_thr, T_thr, MPFR_RNDZ);
        // restar log del primo más grande de la base como margen
        mpfr_sub(T_thr, T_thr, qs_data->base.primes[qs_data->base.length - 1].log_value, MPFR_RNDZ);

        mpz_clears(maxQ, tmp, NULL);
    }
    float T_val = mpfr_get_flt(T_thr, MPFR_RNDZ);
    mpfr_clear(T_thr);

    // Para cada primo p de la base, calcular las raíces de Q(x) ≡ 0 (mod p)
    // Q(x) = a*x² + 2*b*x + c
    // Raíces: x ≡ a⁻¹ * (-b ± sqrt(N)) (mod p)
    // ya que a*x² + 2*b*x + c ≡ 0 (mod p) y b² - a*c = N
    for (long i = 0; i < qs_data->base.length; i++) {
        unsigned long p_ul = mpz_get_ui(qs_data->base.primes[i].value);
        if (p_ul < 2) continue;
        float logp = mpfr_get_flt(qs_data->base.primes[i].log_value, MPFR_RNDZ);

        mpz_t p_mp, r1, r2;
        mpz_inits(p_mp, r1, r2, NULL);
        mpz_set(p_mp, qs_data->base.primes[i].value);

        // Shanks-Tonelli: r1, r2 son raíces de N mod p
        shanksTonelli(qs_data->n, p_mp, r1, r2);

        // Convertir a raíces de Q(x) ≡ 0 (mod p):
        // x ≡ a⁻¹ * (r - b) (mod p)
        mpz_t a_inv, b_mod;
        mpz_inits(a_inv, b_mod, NULL);

        // Si a no es invertible mod p (p divide a), manejar caso especial
        if (mpz_invert(a_inv, qs_data->poly.a, p_mp) == 0) {
            // p | a → Q(x) = 2*b*x + c (mod p), una sola raíz
            // x ≡ -(2b)⁻¹ * c (mod p)
            mpz_t twob;
            mpz_init(twob);
            mpz_mul_ui(twob, qs_data->poly.b, 2);
            mpz_mod(twob, twob, p_mp);
            if (mpz_invert(a_inv, twob, p_mp) != 0) {
                mpz_t s;
                mpz_init(s);
                mpz_mod(s, qs_data->poly.c, p_mp);
                mpz_neg(s, s);
                mpz_mul(s, s, a_inv);
                mpz_mod(s, s, p_mp);
                long start = mpz_get_si(s);
                // mapear a array: pos en S = x + xmax
                // x puede ser start, start+p, start+2p, ...
                // También x negativo: start - p, start - 2p, ...
                for (long x = start; x <= (long)xmax; x += (long)p_ul) {
                    long idx = x + (long)xmax;
                    if (idx >= 0 && idx < (long)total)
                        S[idx] += logp;
                }
                for (long x = start - (long)p_ul; x >= -(long)xmax; x -= (long)p_ul) {
                    long idx = x + (long)xmax;
                    if (idx >= 0 && idx < (long)total)
                        S[idx] += logp;
                }
                mpz_clear(s);
            }
            mpz_clear(twob);
            mpz_clears(p_mp, r1, r2, a_inv, b_mod, NULL);
            continue;
        }

        mpz_mod(b_mod, qs_data->poly.b, p_mp);

        // s1 = a_inv * (r1 - b) mod p
        mpz_t s1, s2;
        mpz_inits(s1, s2, NULL);
        mpz_sub(s1, r1, b_mod);
        mpz_mul(s1, s1, a_inv);
        mpz_mod(s1, s1, p_mp);

        // s2 = a_inv * (r2 - b) mod p
        mpz_sub(s2, r2, b_mod);
        mpz_mul(s2, s2, a_inv);
        mpz_mod(s2, s2, p_mp);

        long sol1 = mpz_get_si(s1);
        long sol2 = mpz_get_si(s2);

        // Cribar s1: desde sol1 hacia ambos lados
        for (long x = sol1; x <= (long)xmax; x += (long)p_ul) {
            S[x + (long)xmax] += logp;
        }
        for (long x = sol1 - (long)p_ul; x >= -(long)xmax; x -= (long)p_ul) {
            S[x + (long)xmax] += logp;
        }

        // Cribar s2 (si es distinta de s1)
        if (sol1 != sol2) {
            for (long x = sol2; x <= (long)xmax; x += (long)p_ul) {
                S[x + (long)xmax] += logp;
            }
            for (long x = sol2 - (long)p_ul; x >= -(long)xmax; x -= (long)p_ul) {
                S[x + (long)xmax] += logp;
            }
        }

        mpz_clears(p_mp, r1, r2, a_inv, b_mod, s1, s2, NULL);
    }

    // Recolectar candidatos que superan el umbral
    // Pre-asignar con estimación conservadora
    unsigned long capacity = total / 10 + 100;
    long *indices = (long *)malloc(capacity * sizeof(long));
    unsigned long count = 0;

    for (unsigned long i = 0; i < total; i++) {
        if (S[i] >= T_val) {
            if (count >= capacity) {
                capacity *= 2;
                indices = (long *)realloc(indices, capacity * sizeof(long));
            }
            indices[count++] = (long)i - (long)xmax;  // valor x real
        }
    }

    free(S);
    *out_indices = indices;
    *out_count = count;
}

