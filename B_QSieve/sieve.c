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
 * @brief Realiza el tamizado para encontrar números Xi que tienen al menos un factor de la base.
 *
 * @param qs_data Estructura de datos que contiene la información necesaria para el tamizado.
 * @param typeSieving Tipo de tamizado (POSITIVE o NEGATIVE) que determina si se busca en los números positivos o negativos.
 * @return Puntero a un array de números Xi que tienen al menos un factor de la base.
 *         El usuario es responsable de liberar la memoria asignada para este array.
 */
float * sievingNaive(qs_struct * qs_data, enum TypeSieving typeSieving) {
    long intervalo = mpz_get_ui(qs_data->intervalo.length);
    // Inicializar un array para almacenar resultados
    float *S = (float*) malloc(intervalo * sizeof(float));
    memset(S, 0, sizeof(float) * intervalo);

    mpz_t n, raizn;
    mpz_inits(n, raizn, NULL);

    // Inicializar n y raizn con los valores correspondientes
    mpz_set(n, qs_data->n);
    mpz_sqrt(raizn, n);

    // int num_threads = 1;

    // if(CORES > 1){
    //     if (typeSieving == POSITIVE) {
    //         // Redondeo hacia arriba
    //         num_threads = (int) ceil(CORES / 2.0);
    //     } else {
    //         // Redondeo hacia abajo
    //         num_threads = (int) floor(CORES / 2.0);
    //     }
    // }

    // printf("Num_T %d\n",num_threads);
    // fflush(stdout);

    #pragma omp parallel for schedule(dynamic) num_threads(CORES)
    for (int i = 0; i < qs_data->base.length; i++) {
        mpz_t x1, x2, p;
        mpz_inits(x1, x2, p, NULL);

        mpz_set(p, qs_data->base.primes[i].value);

        //printf("Iteracion: %d con primo:%ld desde hilo:%d, positivo:%d\n",i, mpz_get_ui(p), omp_get_thread_num(),typeSieving);

        float logp = mpfr_get_flt(qs_data->base.primes[i].log_value, MPFR_RNDZ);

        // Calcular raíces usando el método de Shanks-Tonelli
        shanksTonelli(n, p, x1, x2);
        mpz_sub(x1, x1, raizn);
        mpz_mod(x1, x1, p);

        mpz_sub(x2, x2, raizn);
        mpz_mod(x2, x2, p);

        // Intercambiar x1 y x2 si es necesario
        if (mpz_cmp(x1, x2) == 1) {
            mpz_swap(x1,x2);
        }

        mpz_t d1, d2, Mp;
        mpz_inits(d1, d2, Mp, NULL);

        // Calcular diferencias y límite superior Mp
        mpz_sub(x2, x2, x1);
        mpz_set(d1, x2);
        mpz_sub(d2, p, d1);
        mpz_sub(Mp, qs_data->intervalo.length, d1);

        mpz_t x;
        mpz_init(x);

        // Ajustar x1 para el tamizado NEGATIVE
        if (typeSieving == NEGATIVE) {
            mpz_sub(x1, x1, d2);
            mpz_mul_si(x1,x1,-1);
        }

        // Realizar el tamizado
        for (mpz_set(x, x1); mpz_cmp(x, Mp) == -1;) {
            #pragma omp atomic
                S[mpz_get_ui(x)] = S[mpz_get_ui(x)] + logp;
                mpz_add(x, x, d1);

            #pragma omp atomic
                S[mpz_get_ui(x)] = S[mpz_get_ui(x)] + logp;
                mpz_add(x, x, d2);
        }

        // Ajustar el último elemento si es necesario
        if (mpz_cmp(x, qs_data->intervalo.length) == -1) {
            #pragma omp atomic
                S[mpz_get_ui(x)] = S[mpz_get_ui(x)] + logp;
        }

        mpz_clears(x1, x2, x, d1, d2, Mp, p, NULL);
    }

    mpz_clears(n, raizn, NULL);
    
    return S;
}

/**
 * @brief Realiza el proceso de tamizado para identificar números B-suaves en el intervalo.
 *
 * Algoritmo tomado del articulo Block Sieving Algorithms Georg (Wambach and Hannes Wettig)
 * 
 * @param qs_data Estructura de datos que contiene la información necesaria para el proceso de tamizado.
 * @param length Puntero para almacenar la longitud del array Xi resultante.
 * @return Puntero al array Xi que contiene los números B-suaves encontrados en el intervalo.
 */
mpz_t *sieving(qs_struct *qs_data, unsigned long *length) {

    unsigned long long intervalLength = mpz_get_ui(qs_data->intervalo.length);
    mpz_t *Xi = (mpz_t *)malloc((intervalLength*0.2) * sizeof(mpz_t));
    long contXi = 0;

    mpz_t raizn;
    mpz_init(raizn);
    mpz_sqrt(raizn, qs_data->n);
    // Calcular T = log(sqrt(2N)M)-Delta (Delta = log(ultimo primo base))
    mpfr_t T;
    mpfr_init_set_ui(T, 2, MPFR_RNDZ);
    mpfr_mul_z(T, T, qs_data->n, MPFR_RNDZ);
    mpfr_sqrt(T, T, MPFR_RNDZ);
    mpfr_mul_z(T, T, qs_data->intervalo.length, MPFR_RNDZ);
    mpfr_log(T, T, MPFR_RNDZ);
    mpfr_sub(T, T, qs_data->base.primes[qs_data->base.length - 1].log_value, MPFR_RNDZ);
    unsigned long uT = mpfr_get_ui(T, MPFR_RNDZ);
    float *sp, *sn;

    printf("T:%ld\n",uT);
    
    //Cribado positivo
    unsigned long sizeLocalInterval = (intervalLength*0.1)/CORES;
    double start_time = omp_get_wtime(); // Tiempo de inicio
    sp = sievingNaive(qs_data, POSITIVE);

    #pragma omp parallel num_threads(CORES)
    {
        unsigned long local_contXi = 0;
        mpz_t *local_Xi = (mpz_t *)malloc((sizeLocalInterval) * sizeof(mpz_t));

        #pragma omp for schedule(static)
        for (unsigned long i = 0; i < intervalLength; i++) {
            if (sp[i] > uT) {
                mpz_init(local_Xi[local_contXi]);
                mpz_add_ui(local_Xi[local_contXi],raizn,i);
                local_contXi++;
            }
        }

        #pragma omp critical
        {
            memcpy(&Xi[contXi], local_Xi, local_contXi * sizeof(mpz_t));
            contXi += local_contXi;
        }

        free(local_Xi);
    }

    free(sp);
    double end_time = omp_get_wtime(); // Tiempo de fin
    printf("Tiempo de ejecución de la sección positiva: %f segundos\n", end_time - start_time);

    //Cribado negativo
    start_time = omp_get_wtime(); // Tiempo de inicio
    sn = sievingNaive(qs_data, NEGATIVE);

    #pragma omp parallel num_threads(CORES)
    {
        unsigned long local_contXi = 0;
        mpz_t *local_Xi = (mpz_t *)malloc((sizeLocalInterval) * sizeof(mpz_t));

        #pragma omp for schedule(static)
        for (unsigned long i = 0; i < intervalLength; i++) {
            if (sn[i] > uT) {
                mpz_init(local_Xi[local_contXi]);
                mpz_sub_ui(local_Xi[local_contXi],raizn,i);
                local_contXi++;
            }
        }

        #pragma omp critical
        {
            memcpy(&Xi[contXi], local_Xi, local_contXi * sizeof(mpz_t));
            contXi += local_contXi;
        }

        free(local_Xi);
    }

    free(sn);
    end_time = omp_get_wtime(); // Tiempo de fin
    printf("Tiempo de ejecución de la sección negativa: %f segundos\n", end_time - start_time);

    *length = contXi;

    // Liberar memoria utilizada en el tamizado
    mpz_clear(raizn);
    mpfr_clear(T);

    // Devolver el array Xi resultante
    return Xi;
}

	
