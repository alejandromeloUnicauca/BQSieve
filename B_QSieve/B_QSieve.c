#include <ctype.h>
#include <gmp.h>
#include <getopt.h>
#include <math.h>
#include <mpfr.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structsqs.h"
#include "sieve.h"
#include "polynomial.h"
#include "factoring.h"
#include <time.h>


void createBlocks(int n, qs_struct * qs_data);
void crearMatrizNula(qs_struct * qs_data);
void imprimirMatriz(matrix matriz);  
void getSieveParams(mpz_t n, sieve_param_t *params);
long generatePrimesBase(mpz_t n, long bound, prime * primes);
void freeStruct(qs_struct * qs_data);
void parseArgs(int argc, char **argv, int *flagd, int *flagh, char **hdvalue, char **bvalue, char **cvalue);
void usage();
extern int VERBOSE;

//Cantidad de cores que se quieren usar para el cribado, 1 por defecto
int CORES = 1;
int SIEVE_MULT = 1; /* multiplicador del sieve_size de la tabla msieve */
// Modo verbose: 0 = silencioso (solo resultado + tiempo), 1 = log detallado
int VERBOSE = 0;

/*--------------------------------------------------------------------
 * choose_multiplier — Multiplicador de Knuth-Schroeppel modificado
 *
 * Elige k (squarefree, pequeño) tal que k*N tenga la mayor cantidad
 * de residuos cuadráticos entre primos pequeños → más B-smooth.
 * Algoritmo adaptado de msieve v1.46 (Jason Papadopoulos, dominio público).
 *--------------------------------------------------------------------*/
#define NUM_TEST_PRIMES_KS 300

static const unsigned int ks_mult_list[] = {
    1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19,
    21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38,
    39, 41, 42, 43, 46, 47, 51, 53, 55, 57, 58, 59,
    61, 62, 65, 66, 67, 69, 70, 71, 73
};
#define NUM_KS_MULTS (sizeof(ks_mult_list)/sizeof(ks_mult_list[0]))

static unsigned int choose_multiplier(mpz_t n, unsigned int fb_size) {
    unsigned int i, j;
    unsigned int num_primes;
    double best_score;
    unsigned int best_mult;
    double scores[NUM_KS_MULTS];
    unsigned int num_multipliers;

    /* Usar min(2*fb_size, NUM_TEST_PRIMES_KS) primos para puntuar */
    num_primes = 2 * fb_size;
    if (num_primes > NUM_TEST_PRIMES_KS)
        num_primes = NUM_TEST_PRIMES_KS;

    /* Leer primos del archivo primes.txt */
    FILE *fp = fopen("primes.txt", "r");
    if (!fp) {
        fprintf(stderr, "choose_multiplier: falta primes.txt\n");
        return 1;
    }

    /* Leer hasta num_primes primos en un buffer */
    unsigned long *test_primes = (unsigned long *)malloc(num_primes * sizeof(unsigned long));
    unsigned int n_read = 0;
    char buf[BUFSIZ];
    while (n_read < num_primes && fgets(buf, BUFSIZ, fp) != NULL) {
        char *ptr = buf;
        while (*ptr && (*ptr == ' ' || *ptr == '\t' || *ptr == '\n')) ptr++;
        if (*ptr == '\0') continue;
        unsigned long p = strtoul(ptr, NULL, 10);
        if (p >= 2) {
            test_primes[n_read++] = p;
        }
    }
    fclose(fp);
    num_primes = n_read;

    /* Paso 1: evaluar la contribución del primo 2 y penalizar por tamaño del multiplicador */
    unsigned long n_mod_8 = mpz_fdiv_ui(n, 8);

    double ln2 = log(2.0);

    for (i = 0; i < NUM_KS_MULTS; i++) {
        unsigned int curr_mult = ks_mult_list[i];
        unsigned int knmod8 = (unsigned int)((curr_mult * n_mod_8) % 8);
        double logmult = log((double)curr_mult);

        /* Penalización: multiplicadores grandes hacen k*N más grande */
        scores[i] = 0.5 * logmult;

        /* Bonus por el primo 2 según k*N mod 8 */
        switch (knmod8) {
        case 1: scores[i] -= 2 * ln2; break;
        case 5: scores[i] -= ln2; break;
        case 3:
        case 7: scores[i] -= 0.5 * ln2; break;
        /* knmod8 par: no hay bonus (multiplicadores pares empiezan con desventaja) */
        }
    }
    num_multipliers = NUM_KS_MULTS;

    /* Paso 2: para cada primo p de test, evaluar contribución log(p)/(p-1) */
    for (i = 1; i < num_primes; i++) {  /* empezar en 1 para saltar p=2 */
        unsigned long prime = test_primes[i];
        double contrib = log((double)prime) / (double)(prime - 1);
        unsigned long n_mod_p = mpz_fdiv_ui(n, prime);

        for (j = 0; j < num_multipliers; j++) {
            unsigned int curr_mult = ks_mult_list[j];
            unsigned long kn_mod_p = (n_mod_p * (curr_mult % prime)) % prime;

            /* Si k*N es residuo cuadrático mod prime (o prime | k*N) */
            if (kn_mod_p == 0) {
                /* prime divide k*N → solo una raíz */
                scores[j] -= contrib;
            } else {
                /* Legendre symbol: kn_mod_p^((p-1)/2) mod p */
                mpz_t tmp_kn, tmp_p;
                mpz_inits(tmp_kn, tmp_p, NULL);
                mpz_set_ui(tmp_kn, kn_mod_p);
                mpz_set_ui(tmp_p, prime);
                int leg = mpz_legendre(tmp_kn, tmp_p);
                mpz_clears(tmp_kn, tmp_p, NULL);

                if (leg == 1) {
                    /* Dos raíces → doble contribución */
                    scores[j] -= 2.0 * contrib;
                }
            }
        }
    }

    free(test_primes);

    /* Paso 3: elegir el multiplicador con mejor score (más negativo = mejor) */
    best_score = 1000.0;
    best_mult = 1;
    for (i = 0; i < num_multipliers; i++) {
        if (scores[i] < best_score) {
            best_score = scores[i];
            best_mult = ks_mult_list[i];
        }
    }
    return best_mult;
}

int main(int argc, char **argv)
{
	int flagd = 0; 
	int flagh = 0;
	char *hdvalue = NULL;
	char *bvalue = NULL;
	char *cvalue = NULL;
	
	parseArgs(argc, argv, &flagd, &flagh, &hdvalue, &bvalue, &cvalue);
	
	//Declaracion de variables
	qs_struct qs_data;
	clock_t t_inicio, t_final;

	//Inicializar campos a valores seguros
	qs_data.n_BSuaves = 0;
	qs_data.base.primes = NULL;
	qs_data.base.sp = NULL;
	qs_data.base.length = 0;
	qs_data.blocks.block = NULL;
	qs_data.blocks.length = 0;
	qs_data.mat.data = NULL;
	qs_data.mat.n_rows = 0;
	qs_data.mat.n_cols = 0;
	qs_data.intervalo.Xi = NULL;
	qs_data.intervalo.Qxi = NULL;
	qs_data.intervalo.length_Xi = 0;
	qs_data.intervalo.length_Qxi = 0;
	// inicializar mpz_t del polinomio
	mpz_init(qs_data.poly.a);
	mpz_init(qs_data.poly.b);
	mpz_init(qs_data.poly.c);
	// inicializar roota persistente para MPQS (0 indica no inicializado)
	mpz_init(qs_data.roota);
	mpz_set_ui(qs_data.roota, 0);
	// inicializar estado SIQS
	memset(&qs_data.siqs_state, 0, sizeof(siqs_poly_state));
	// inicializar tabla de parciales
	qs_data.partials.entries = NULL;
	qs_data.partials.n = 0;
	qs_data.partials.capacity = 0;
	qs_data.large_prime_bound = 0;
	qs_data.large_prime_bound2 = 0;
	qs_data.max_fb2 = 0;
	qs_data.n_dlp_stored = 0;
	qs_data.n_dlp_combined = 0;
	qs_data.multiplier = 1;
	mpz_inits(qs_data.n,qs_data.intervalo.length,NULL);
	if(bvalue!=NULL)qs_data.blocks.length = atol(bvalue);
	else qs_data.blocks.length = 0;
	
	if(cvalue!=NULL)
		CORES = atoi(cvalue);
	
	//si se usa el flag -d se asigna un numero decimal 
	//si se usa el flag -h se asigna un numero hexadecimal
	//si no se termina
	if(flagd == 1){
		if(mpz_set_str(qs_data.n, hdvalue, 10)==-1){
			fprintf(stderr,"N no es un numero valido\n");
			usage();
			exit(EXIT_FAILURE);
		}
	}else if(flagh == 1){
		if(mpz_set_str(qs_data.n, hdvalue, 16)==-1){
			fprintf(stderr,"N no es un numero valido\n");
			usage();
			exit(EXIT_FAILURE);
		}
	}else{
		usage();
		exit(EXIT_FAILURE);
	}
	
	if(VERBOSE) gmp_printf("N:%Zd\n", qs_data.n);
	
	//digitos de N
	size_t sizeN = mpz_sizeinbase(qs_data.n, 10);
	if(VERBOSE) printf("Numero de digitos decimales: %zu\n",sizeN);

	if(VERBOSE) printf("Cores:%d\n",CORES);
	
	//Obtener parámetros de criba interpolados según el tamaño de N (antes de multiplicar por k)
	getSieveParams(qs_data.n, &qs_data.sieve_params);
	qs_data.base.length = qs_data.sieve_params.fb_size;
	if(VERBOSE) printf("Parámetros de criba: bits=%u, fb_size=%u, sieve_size=%u, large_mult=%u\n",
		qs_data.sieve_params.bits, qs_data.sieve_params.fb_size,
		qs_data.sieve_params.sieve_size, qs_data.sieve_params.large_mult);

	// Elegir multiplicador Knuth-Schroeppel
	qs_data.multiplier = choose_multiplier(qs_data.n, qs_data.sieve_params.fb_size);
	if(VERBOSE) printf("Multiplicador Knuth-Schroeppel: k=%u\n", qs_data.multiplier);
	if (qs_data.multiplier > 1) {
		mpz_mul_ui(qs_data.n, qs_data.n, qs_data.multiplier);
		if(VERBOSE) gmp_printf("kN:%Zd\n", qs_data.n);
	}

	if(VERBOSE) printf("Longitud de la base de primos:%ld\n", qs_data.base.length);
	
	//Generar base de primos
	qs_data.base.primes = (prime*)malloc((qs_data.base.length)*sizeof(prime));
	
	t_inicio = clock();
	if(VERBOSE) printf("Generando base de primos...\n");
	long residuos = generatePrimesBase(qs_data.n,qs_data.base.length,qs_data.base.primes);
	t_final = clock();
	// ajustar longitud real de la base al número de residuos encontrados
	qs_data.base.length = residuos;
	// reducir buffer al tamaño real
	if (residuos > 0) {
		qs_data.base.primes = (prime*)realloc(qs_data.base.primes, residuos * sizeof(prime));
	}
	if(VERBOSE) printf("Base de primos generada. %ld primos en la base\n",residuos);

	// Precomputar raíces sqrt(N) mod p y campos nativos (uint32/uint8)
	sieve_precompute_roots(&qs_data);

	//Intervalo de criba: sieve_size de la tabla msieve × SIEVE_MULT (flag -s)
	mpz_set_ui(qs_data.intervalo.length,
	           (unsigned long)qs_data.sieve_params.sieve_size * (unsigned long)SIEVE_MULT);

	// Calcular large_prime_bound = large_mult * primo_más_grande_de_la_base
	{
		unsigned long p_max = mpz_get_ui(qs_data.base.primes[qs_data.base.length - 1].value);
		qs_data.large_prime_bound = (unsigned long)qs_data.sieve_params.large_mult * p_max;
		if(VERBOSE) printf("Large prime bound: %lu (large_mult=%u × p_max=%lu)\n",
			qs_data.large_prime_bound, qs_data.sieve_params.large_mult, p_max);

		/* max_fb2 = p_max² (umbral mínimo: cofactores < max_fb2 son primos → 1LP) */
		qs_data.max_fb2 = (unsigned long long)p_max * (unsigned long long)p_max;

		/* 2LP habilitado para N >= 85 dígitos (~282 bits), como msieve.
		 * large_prime_bound2 = LP_bound^1.8 (cofactor máximo para 2LP).
		 * Ambos factores del cofactor deben ser < LP_bound. */
		unsigned int nbits = (unsigned int)mpz_sizeinbase(qs_data.n, 2);
		if (nbits >= 282 && qs_data.base.length >= 800) {
			double lp = (double)qs_data.large_prime_bound;
			qs_data.large_prime_bound2 = (unsigned long long)(lp * pow(lp, 0.8));
			if(VERBOSE) printf("Double Large Prime bound: %llu (%u-%llu bits)\n",
				qs_data.large_prime_bound2,
				(unsigned int)(log2((double)qs_data.max_fb2)),
				(unsigned long long)(log2((double)qs_data.large_prime_bound2)));
		} else {
			qs_data.large_prime_bound2 = 0;
			if(VERBOSE) printf("Double Large Primes: deshabilitado (N < 282 bits o fb < 800)\n");
		}
	}

	double segundos = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
	if(VERBOSE) printf("tiempo de creacion de la base:%fs\n",segundos);
	
	//Crear bloques de la base
	if(qs_data.blocks.length > 0){
		int blockLength = ceil((float)qs_data.base.length/qs_data.blocks.length);
		if(VERBOSE) printf("Creando bloques...\n");
		t_inicio = clock();
		createBlocks(blockLength,&qs_data);
		t_final = clock();
		if(VERBOSE) printf("Bloques creados: %ld\n",qs_data.blocks.length);
		double segundos = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
		if(VERBOSE) printf("tiempo de creacion de los bloques:%fs\n",segundos);	
	}
	
	if(VERBOSE) printf("Cribando...\n");
	double start_time = omp_get_wtime();

    // xmax define el rango de criba [-xmax..+xmax], tomado del intervalo del polinomio
    unsigned long xmax = mpz_get_ui(qs_data.intervalo.length);

    double end_time = omp_get_wtime();
    double segundosCriba = end_time - start_time;
    if(VERBOSE) printf("xmax: %lu\n", xmax);
    if(VERBOSE) printf("tiempo de preparación de criba: %f segundos\n", segundosCriba);

	if(VERBOSE) printf("Calculando Polinomio...\n");
	t_inicio = clock();
	double t_poly_wall0 = omp_get_wtime();
	crearMatrizNula(&qs_data);
	polinomio_open();

	int res = 1;
	qs_data.intervalo.Qxi = NULL;

	long polinomio_count = 0;
	long prev_n_BSuaves = qs_data.n_BSuaves;
	long full_relations = 0;

	double t_gen = 0, t_sieve = 0, t_evalQ = 0, t_factor = 0;
	long n_cand_total = 0;

	while(res==1){
		double _t0, _t1;

		// MPQS: generar nuevo polinomio
		_t0 = omp_get_wtime();
		generate_mpqs_poly(&qs_data);
		_t1 = omp_get_wtime();
		t_gen += _t1 - _t0;

		// === CRIBA LOGARÍTMICA ===
		// Pre-filtrar candidatos con criba logarítmica antes de trial division
		long *sieve_candidates = NULL;
		unsigned long n_candidates = 0;
		_t0 = omp_get_wtime();
		sieve_mpqs(&qs_data, xmax, &sieve_candidates, &n_candidates);
		_t1 = omp_get_wtime();
		t_sieve += _t1 - _t0;

		if (n_candidates == 0) {
			free(sieve_candidates);
			polinomio_count++;
			continue;
		}

		// Construir Xi y Qxi solo para los candidatos de la criba
		unsigned long npos = n_candidates;
		n_cand_total += (long)npos;

		// Liberar Xi y Qxi previos
		if (qs_data.intervalo.Xi != NULL) {
			for (unsigned long i = 0; i < qs_data.intervalo.length_Xi; i++)
				mpz_clear(qs_data.intervalo.Xi[i]);
			free(qs_data.intervalo.Xi);
			qs_data.intervalo.Xi = NULL;
		}
		if (qs_data.intervalo.Qxi != NULL) {
			for (unsigned long i = 0; i < qs_data.intervalo.length_Qxi; i++)
				mpz_clear(qs_data.intervalo.Qxi[i]);
			free(qs_data.intervalo.Qxi);
			qs_data.intervalo.Qxi = NULL;
		}

		// Asignar nuevos arrays
		qs_data.intervalo.Xi = (mpz_t *)malloc(npos * sizeof(mpz_t));
		qs_data.intervalo.Qxi = (mpz_t *)malloc(npos * sizeof(mpz_t));
		qs_data.intervalo.length_Xi = npos;
		qs_data.intervalo.length_Qxi = npos;

		_t0 = omp_get_wtime();
		for (unsigned long i = 0; i < npos; i++) {
			mpz_init(qs_data.intervalo.Xi[i]);
			mpz_set_si(qs_data.intervalo.Xi[i], sieve_candidates[i]);
			mpz_init(qs_data.intervalo.Qxi[i]);
			eval_mpqs_Qx(&qs_data, qs_data.intervalo.Xi[i], qs_data.intervalo.Qxi[i]);
		}
		_t1 = omp_get_wtime();
		t_evalQ += _t1 - _t0;

		free(sieve_candidates);

		polinomio_count++;
		_t0 = omp_get_wtime();
		if (qs_data.blocks.length > 0) {
			res = factoringBlocks(&qs_data, npos, 0, xmax);
		} else {
			res = factoringTrial(&qs_data, npos, 0, xmax);
		}
		_t1 = omp_get_wtime();
		t_factor += _t1 - _t0;

		long found_this = qs_data.n_BSuaves - prev_n_BSuaves;
		if (found_this > 0) {
			/* Clasificar: full vs combined */
			full_relations += found_this; /* ajustado abajo */
			if(VERBOSE) {
				printf("Polinomio %ld: %lu candidatos criba → %ld B_suaves (total: %ld, parciales: %lu)\n",
					polinomio_count, npos, found_this, qs_data.n_BSuaves, qs_data.partials.n);
				fflush(stdout);
			}
		}
		prev_n_BSuaves = qs_data.n_BSuaves;
	}
	polinomio_close();
	if(VERBOSE) printf("Polinomios procesados: %ld\n", polinomio_count);
	if(VERBOSE) fflush(stdout);
	if(VERBOSE) printf("Numeros B_Suaves encontrados:%ld\n",qs_data.n_BSuaves);
	if(VERBOSE) printf("Parciales almacenadas sin emparejar: %lu\n", qs_data.partials.n);
	if (VERBOSE && qs_data.large_prime_bound2 > 0) {
		printf("2LP: almacenadas=%lu, combinadas=%lu\n",
			qs_data.n_dlp_stored, qs_data.n_dlp_combined);
	}
	t_final = clock();
	double t_poly_wall = omp_get_wtime() - t_poly_wall0;
	double segundosPolinomio = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
	if(VERBOSE) printf("tiempo de calculo del polinomio:%fs (CPU) | %fs (wall)\n",
	                   segundosPolinomio, t_poly_wall);
	if (VERBOSE && t_poly_wall > 0) {
		double tw = t_poly_wall;
		printf("  ├─ generate_poly: %.3fs (%5.1f%%)\n", t_gen,    100.0*t_gen/tw);
		printf("  ├─ sieve_mpqs:    %.3fs (%5.1f%%)\n", t_sieve,  100.0*t_sieve/tw);
		printf("  ├─ eval_Qx:       %.3fs (%5.1f%%) — %ld candidatos totales\n",
		       t_evalQ, 100.0*t_evalQ/tw, n_cand_total);
		printf("  └─ trial+combine: %.3fs (%5.1f%%)\n", t_factor, 100.0*t_factor/tw);
		double sum = t_gen + t_sieve + t_evalQ + t_factor;
		printf("  (suma sub-fases: %.3fs = %.1f%% del wall; overhead/otros: %.3fs)\n",
		       sum, 100.0*sum/tw, tw - sum);
		extern void print_factoring_stats(double total_wall);
		print_factoring_stats(t_factor);
		extern void print_sieve_stats(double total_wall);
		print_sieve_stats(t_sieve);
	}
	
	
	if(VERBOSE) printf("Escribiendo matriz...");
	imprimirMatriz(qs_data.mat);

	// Guardar roota en archivo para que mulPoli lo use
	FILE *fr = fopen("roota.txt", "w");
	if (fr) {
		mpz_out_str(fr, 10, qs_data.roota);
		fprintf(fr, "\n");
		fclose(fr);
	}

	// Guardar multiplicador para que mulPoli divida los factores espúreos
	if (qs_data.multiplier > 1) {
		FILE *fm = fopen("multiplier.txt", "w");
		if (fm) {
			fprintf(fm, "%u\n", qs_data.multiplier);
			fclose(fm);
		}
	}
	
	//Liberar Memoria
	freeStruct(&qs_data); 
	
	exit(EXIT_SUCCESS);
}

/**
 * @brief Imprime una matriz dispersa en un archivo.
 *
 * Esta función toma una matriz dispersa representada por la estructura `matrix` y la imprime
 * en un archivo llamado "matrix.txt". La matriz dispersa se representa indicando el número de
 * elementos no cero en cada fila, seguido de los índices de las columnas donde se encuentran
 * los elementos no cero.
 *
 * @param matriz Estructura que representa la matriz dispersa.
 */
void imprimirMatriz(matrix matriz) {
    FILE* f = fopen("matrix.txt", "w");  // Abrir el archivo en modo escritura

    // Imprimir las dimensiones de la matriz en la primera línea del archivo
    fprintf(f, "%d %d\n", matriz.n_rows, matriz.n_cols);

    int v[matriz.n_cols];  // Arreglo para almacenar índices de elementos no cero
    int cont;              // Contador de elementos no cero

    // Recorrer filas de la matriz
    for (int i = 0; i < matriz.n_rows; i++) {
        cont = 0;

        // Recorrer columnas de la matriz
        for (int j = 0; j < matriz.n_cols; j++) {
            // Si el elemento en la posición (i, j) es 1, almacenar el índice j
            if (matriz.data[i][j] == 1) {
                v[cont] = j;
                cont++;
            }
        }

        // Imprimir el número de elementos no cero en la fila y sus índices
        fprintf(f, "%d ", cont);
        for (int k = 0; k < cont; k++) {
            fprintf(f, "%d ", v[k]);
        }

        fprintf(f, "\n");  // Nueva línea para la siguiente fila
    }

    fclose(f);  // Cerrar el archivo
}

void crearMatrizNula(qs_struct * qs_data){
	// Necesitamos más filas (relaciones) que columnas para que el espacio nulo
	// tenga dimensión suficiente. extra relaciones dan más soluciones independientes.
	int extra = 64; // al menos 64 relaciones extra para tener ~64 soluciones
	qs_data->mat.n_rows = qs_data->base.length + 1 + extra;
	qs_data->mat.n_cols = qs_data->base.length + 1; // +1 para columna de signo (-1)
	
	//reservar memoria para matriz
	qs_data->mat.data = (int**)malloc(qs_data->mat.n_rows*sizeof(int*));
	   
	for (int i = 0; i < qs_data->mat.n_rows; i++) 
	{
		qs_data->mat.data[i] = (int*)malloc(qs_data->mat.n_cols*sizeof(int));
		memset(qs_data->mat.data[i],0,qs_data->mat.n_cols * sizeof(int));
	}
}

/**
 * @brief crea un archivo bloques.txt donde se almacena la multiplicacion
 * de los residuos separados en bloques de tamaño n
 * @param n:tamaño de los bloques
 */
void createBlocks(int n, qs_struct * qs_data){
	/* Número real de bloques que se necesitan para cubrir toda la base.
	 * El user pidió blocks.length grupos, pero si base.length no es múltiplo
	 * de n el último grupo es más pequeño y los siguientes serían vacíos. */
	unsigned long real_blocks = (qs_data->base.length + n - 1) / n;
	if (real_blocks < qs_data->blocks.length) qs_data->blocks.length = real_blocks;

	qs_data->blocks.block = (prime_block*)malloc(qs_data->blocks.length * sizeof(prime_block));
	for (unsigned long i = 0; i < qs_data->blocks.length; i++)
		qs_data->blocks.block[i].factors = (prime*)malloc(n * sizeof(prime));

	mpz_t mulTemp;
	mpz_init(mulTemp);
	mpz_set_ui(mulTemp, 1);

	int contBlock = 0;
	int contFact = 0;
	for (int i = 0; i < qs_data->base.length; i++) {
		if (contFact == n) {
			qs_data->blocks.block[contBlock].length = contFact;
			mpz_init(qs_data->blocks.block[contBlock].prod_factors);
			mpz_set(qs_data->blocks.block[contBlock].prod_factors, mulTemp);
			mpz_set_ui(mulTemp, 1);
			contBlock++;
			contFact = 0;
		}
		mpz_init(qs_data->blocks.block[contBlock].factors[contFact].value);
		mpz_set(qs_data->blocks.block[contBlock].factors[contFact].value, qs_data->base.primes[i].value);
		mpz_mul(mulTemp, mulTemp, qs_data->base.primes[i].value);
		contFact++;
	}
	qs_data->blocks.block[contBlock].length = contFact;
	mpz_init(qs_data->blocks.block[contBlock].prod_factors);
	mpz_set(qs_data->blocks.block[contBlock].prod_factors, mulTemp);
	mpz_clear(mulTemp);
}

void freeStruct(qs_struct * qs_data){
	
	//liberar memoria de la base
	for (int i = 0; i < qs_data->base.length; i++)
	{
		//gmp_printf("P:%Zd",qs_data->base.primes[i].value);
		mpz_clear(qs_data->base.primes[i].value);
		//mpfr_printf ("log(p):%.2Rf\n", qs_data->base.primes[i].log_value);
		mpfr_clear(qs_data->base.primes[i].log_value);
	}

	free(qs_data->base.primes);
	free(qs_data->base.sp);

	//liberar memoria de los bloques
	if(qs_data->blocks.length > 0){
		for (int i = 0; i < qs_data->blocks.length ; i++)
		{
			//printf("Bloque %d:",i);
			for (int j = 0; j < qs_data->blocks.block[i].length; j++)
			{
				//gmp_printf("%Zd,",qs_data->blocks.block[i].factors[j].value);
				mpz_clear(qs_data->blocks.block[i].factors[j].value);
			}
			//gmp_printf("%Zd,",qs_data->blocks.block[i].prod_factors);
			//printf("\n");
			mpz_clear(qs_data->blocks.block[i].prod_factors);
			free(qs_data->blocks.block[i].factors);
		}
		free(qs_data->blocks.block);
	}

	if(qs_data->intervalo.Qxi!=NULL){
		unsigned long long lenQ = qs_data->intervalo.length_Qxi;
		for (unsigned long long i = 0; i < lenQ; i++)
		{
			mpz_clear(qs_data->intervalo.Qxi[i]);
		}
		free(qs_data->intervalo.Qxi);
		qs_data->intervalo.Qxi = NULL;
		qs_data->intervalo.length_Qxi = 0;
	}

	if(qs_data->intervalo.Xi!=NULL){
		unsigned long long lenXi = qs_data->intervalo.length_Xi;
		for (unsigned long long i = 0; i < lenXi; i++)
		{
			mpz_clear(qs_data->intervalo.Xi[i]);
		}
		free(qs_data->intervalo.Xi);
		qs_data->intervalo.Xi = NULL;
		qs_data->intervalo.length_Xi = 0;
	}

	// liberar parciales
	for (unsigned long long i = 0; i < qs_data->partials.n; i++) {
		mpz_clear(qs_data->partials.entries[i].lhs);
		mpz_clear(qs_data->partials.entries[i].Qx);
		mpz_clear(qs_data->partials.entries[i].roota);
		mpz_clear(qs_data->partials.entries[i].a_value);
		free(qs_data->partials.entries[i].exponents);
	}
	free(qs_data->partials.entries);

	// liberar roota y polinomio MPQS
	mpz_clear(qs_data->roota);
	mpz_clear(qs_data->poly.a);
	mpz_clear(qs_data->poly.b);
	mpz_clear(qs_data->poly.c);

	// liberar estado SIQS si fue inicializado
	if (qs_data->siqs_state.initialized) {
		for (unsigned int i = 0; i < MAX_SIQS_FACTORS; i++) {
			mpz_clear(qs_data->siqs_state.factors[i]);
			mpz_clear(qs_data->siqs_state.Bvals[i]);
		}
		mpz_clear(qs_data->siqs_state.target_a);
	}

	
	//liberar memoria de la matriz
	   
	for (int i = 0; i < qs_data->mat.n_rows; i++) 
	{
		free(qs_data->mat.data[i]);
	}
	
	free(qs_data->mat.data);
}


/**
 * @brief calcula los residuos cuadraticos del numero n usando la funcion legendre,
 * los primos se obtienen del archivo primes.txt y los almacena en un archivo residuos.txt 
 * @param n: numero que se le evaluaran los residuos cuadraticos
 * @param base_length: numero de residuos que se requieren para el numero n
 * @return retorna el numero de 
 * residuos encontrados
 */
long generatePrimesBase(mpz_t n, long bound, prime * primes){
    long contRes = 0; // contador de residuos encontrados
    long contPrimos = 0; // contador de primos leídos del archivo

    mpz_t p; // variable temporal para los primos del archivo
    mpz_init(p);

    FILE * file; // file primes
    // si el archivo primes.txt no existe termina
    if ((file = fopen("primes.txt", "r")) == NULL) // open file
    {
        fprintf(stderr,"Falta archivo primes.txt\n");
        exit(EXIT_FAILURE);
    }

    char buf[BUFSIZ];
    while (fgets(buf, BUFSIZ, file) != NULL) {
        // eliminar espacios en blanco iniciales
        char *ptr = buf;
        while (*ptr && isspace((unsigned char)*ptr)) ptr++;
        if (*ptr == '\0') continue;

        if (mpz_set_str(p, ptr, 10) != 0) continue; // parse error

        contPrimos++;

        // si n es residuo cuadratico mod p se agrega
        if ((mpz_legendre(n,p) == 1) || (mpz_cmp_ui(p,2) == 0)){
            // asigno memoria a los valores de prime
            mpz_init(primes[contRes].value);
            mpfr_init(primes[contRes].log_value);

            // almaceno el primo y el logaritmo del primo
            mpz_set(primes[contRes].value,p);

            mpfr_t pTemp;
            mpfr_init(pTemp);
            mpfr_set_z(pTemp,p,MPFR_RNDZ);
            mpfr_log(primes[contRes].log_value, pTemp, MPFR_RNDZ); // ln(p)
            primes[contRes].llog_value = mpfr_get_ui(primes[contRes].log_value,MPFR_RNDZ);

            mpfr_clear(pTemp);

            contRes++;

            // parar cuando tengamos suficientes residuos cuadráticos
            if (contRes >= bound) break;
        }
    }

    mpz_clear(p);
    fclose(file);
    if(VERBOSE) gmp_printf("Primo mas grande en la base: %Zd\n", primes[contRes-1].value);
    if(VERBOSE) printf("Primos leidos del archivo: %ld, residuos cuadraticos: %ld de %ld requeridos\n",
           contPrimos, contRes, bound);
    return contRes;
}

/**
 * @brief Tabla de parámetros de criba precompilados, indexados por bits de N.
 * Adaptada de msieve v1.46 (Jason Papadopoulos, dominio público).
 * {bits, fb_size, large_mult, sieve_size}
 */
static const sieve_param_t prebuilt_params[] = {
	{ 64,    100,  40,  1 * 65536},
	{128,    450,  40,  1 * 65536},
	{183,   2000,  40,  1 * 65536},
	{200,   3000,  50,  1 * 65536},
	{212,   5400,  50,  3 * 65536},
	{233,  10000, 100,  3 * 65536},
	{249,  27000, 100,  3 * 65536},
	{266,  50000, 100,  3 * 65536},
	{283,  55000,  80,  3 * 65536},
	{298,  60000,  80,  9 * 65536},
	{315,  80000, 150,  9 * 65536},
	{332, 100000, 150,  9 * 65536},
	{348, 140000, 150,  9 * 65536},
	{363, 210000, 150, 13 * 65536},
	{379, 300000, 150, 17 * 65536},
	{395, 400000, 150, 21 * 65536},
	{415, 500000, 150, 25 * 65536},
	{440, 700000, 150, 33 * 65536},
	{465, 900000, 150, 50 * 65536},
	{490,1100000, 150, 75 * 65536},
	{512,1300000, 150,100 * 65536},
};
#define NUM_PREBUILT_PARAMS (sizeof(prebuilt_params)/sizeof(sieve_param_t))

/**
 * @brief Obtiene los parámetros de criba interpolados según el tamaño en bits de N.
 *
 * Si N cae entre dos entradas de la tabla, se interpola linealmente.
 * Garantiza fb_size >= 100.
 *
 * @param n      Número a factorizar
 * @param params Estructura de salida con los parámetros
 */
void getSieveParams(mpz_t n, sieve_param_t *params) {
	unsigned int bits = (unsigned int)mpz_sizeinbase(n, 2);

	/* Si es más pequeño que la primera entrada, usar la primera */
	if (bits <= prebuilt_params[0].bits) {
		*params = prebuilt_params[0];
		params->bits = bits;
		return;
	}

	/* Si es más grande que la última entrada, usar la última */
	if (bits >= prebuilt_params[NUM_PREBUILT_PARAMS - 1].bits) {
		*params = prebuilt_params[NUM_PREBUILT_PARAMS - 1];
		params->bits = bits;
		return;
	}

	/* Buscar las dos entradas entre las que cae bits */
	unsigned int i;
	for (i = 0; i < NUM_PREBUILT_PARAMS - 1; i++) {
		if (bits < prebuilt_params[i + 1].bits)
			break;
	}

	/* Interpolación lineal ponderada */
	const sieve_param_t *low  = &prebuilt_params[i];
	const sieve_param_t *high = &prebuilt_params[i + 1];
	unsigned int dist = high->bits - low->bits;
	unsigned int wi = bits - low->bits;     /* peso hacia high */
	unsigned int wj = high->bits - bits;    /* peso hacia low  */

	params->bits = bits;
	params->fb_size = (unsigned int)(
		((double)low->fb_size * wj + (double)high->fb_size * wi) / dist + 0.5);
	params->large_mult = (unsigned int)(
		((double)low->large_mult * wj + (double)high->large_mult * wi) / dist + 0.5);
	params->sieve_size = (unsigned int)(
		((double)low->sieve_size * wj + (double)high->sieve_size * wi) / dist + 0.5);

	/* Mínimo de 100 primos en la base */
	if (params->fb_size < 100)
		params->fb_size = 100;
}

void parseArgs(int argc, char **argv, int *flagd, int *flagh, char **hdvalue, char **bvalue, char **cvalue){
	int c;
	opterr = 0;
	
	while ((c = getopt(argc, argv, "d:h:b:c:s:v")) != -1){
		switch(c){
			case 'd':
				if(*flagh == 1){
					fprintf (stderr, "Solo puedes usar -h o -d pero no ambos\n");
					usage();
					exit(EXIT_FAILURE);
				}
				*flagd = 1;
				*hdvalue = optarg;
				break;
			case 'h':
				if(*flagd == 1){
					fprintf (stderr, "Solo puedes usar -h o -d pero no ambos\n");
					usage();
					exit(EXIT_FAILURE);
				}
				*flagh = 1;
				*hdvalue = optarg;
				break;
			case 'b':
				*bvalue = optarg;
				break;
			case 'c':
				*cvalue = optarg;
				break;
			case 's':
				SIEVE_MULT = atoi(optarg);
				if (SIEVE_MULT < 1) SIEVE_MULT = 1;
				break;
			case 'v':
				VERBOSE = 1;
				break;
			case '?':
				if (strchr("h", optopt) != NULL)
					fprintf (stderr, "La opción -%c requiere un argumento.\n", optopt);
				else if (strchr("d", optopt) != NULL)
					fprintf (stderr, "La opción -%c requiere un argumento.\n", optopt);
				else if (strchr("b", optopt) != NULL)
					fprintf (stderr, "La opción -%c requiere un argumento.\n", optopt);
				else if (strchr("c", optopt) != NULL)
					fprintf (stderr, "La opción -%c requiere un argumento.\n", optopt);
				else if (isprint (optopt))
					fprintf(stderr, "Opción desconocida'-%c'.\n", optopt);
				else
					fprintf(stderr, "Carácter no válido '\\x%x'.\n", optopt);
				usage();
				exit(EXIT_FAILURE);
				
				break;
		}
	}
	
	if(argc < 2){
		usage();
		exit(EXIT_FAILURE);
	}
}

void usage(){
	fprintf(stderr,"Uso: ./B_QSieve (-d | -h) <N> [-b <NBLOCKS>] [-c <NCORES>] [-v]\n");
	fprintf(stderr,"Opciones:\n");
	fprintf(stderr,"-d	# Especifica que el numero N es decimal\n");
	fprintf(stderr,"-h	# Especifica que el numero N es hexadecimal\n");
	fprintf(stderr,"-b	# Al usar esta opcion se deben especificar el numero de Bloques\n");
	fprintf(stderr,"-c	# Especifica el numero de procesadores logicos que se quieren usar\n");
	fprintf(stderr,"-v	# Modo verbose: muestra el log detallado de ejecución\n");
}
