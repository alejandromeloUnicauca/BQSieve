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

//Cantidad de cores que se quieren usar para el cribado, 1 por defecto
int CORES = 1;

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
	
	gmp_printf("N:%Zd\n", qs_data.n);
	
	//digitos de N
	size_t sizeN = mpz_sizeinbase(qs_data.n, 10);
	printf("Numero de digitos decimales: %zu\n",sizeN);

	printf("Cores:%d\n",CORES);
	
	//Obtener parámetros de criba interpolados según el tamaño de N
	getSieveParams(qs_data.n, &qs_data.sieve_params);
	qs_data.base.length = qs_data.sieve_params.fb_size;
	printf("Parámetros de criba: bits=%u, fb_size=%u, sieve_size=%u, large_mult=%u\n",
		qs_data.sieve_params.bits, qs_data.sieve_params.fb_size,
		qs_data.sieve_params.sieve_size, qs_data.sieve_params.large_mult);
	printf("Longitud de la base de primos:%ld\n", qs_data.base.length);
	
	//Generar base de primos
	qs_data.base.primes = (prime*)malloc((qs_data.base.length)*sizeof(prime));
	
	t_inicio = clock();
	printf("Generando base de primos...\n");
	long residuos = generatePrimesBase(qs_data.n,qs_data.base.length,qs_data.base.primes);
	t_final = clock();
	// ajustar longitud real de la base al número de residuos encontrados
	qs_data.base.length = residuos;
	// reducir buffer al tamaño real
	if (residuos > 0) {
		qs_data.base.primes = (prime*)realloc(qs_data.base.primes, residuos * sizeof(prime));
	}
	printf("Base de primos generada. %ld primos en la base\n",residuos);

	// Precomputar raíces sqrt(N) mod p y campos nativos (uint32/uint8)
	sieve_precompute_roots(&qs_data);

	//Intervalo de criba: usar sieve_size de la tabla de parámetros
	mpz_set_ui(qs_data.intervalo.length, qs_data.sieve_params.sieve_size);

	// Calcular large_prime_bound = large_mult * primo_más_grande_de_la_base
	{
		unsigned long p_max = mpz_get_ui(qs_data.base.primes[qs_data.base.length - 1].value);
		qs_data.large_prime_bound = (unsigned long)qs_data.sieve_params.large_mult * p_max;
		printf("Large prime bound: %lu (large_mult=%u × p_max=%lu)\n",
			qs_data.large_prime_bound, qs_data.sieve_params.large_mult, p_max);
	}

	double segundos = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
	printf("tiempo de creacion de la base:%fs\n",segundos);
	
	//Crear bloques de la base
	if(qs_data.blocks.length > 0){
		int blockLength = ceil((float)qs_data.base.length/qs_data.blocks.length);
		printf("Creando bloques...\n");
		t_inicio = clock();
		createBlocks(blockLength,&qs_data);
		t_final = clock();
		printf("Bloques creados: %ld\n",qs_data.blocks.length);
		double segundos = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
		printf("tiempo de creacion de los bloques:%fs\n",segundos);	
	}
	
	printf("Cribando...\n");
	double start_time = omp_get_wtime();

    // xmax define el rango de criba [-xmax..+xmax], tomado del intervalo del polinomio
    unsigned long xmax = mpz_get_ui(qs_data.intervalo.length);

    double end_time = omp_get_wtime();
    double segundosCriba = end_time - start_time;
    printf("xmax: %lu\n", xmax);
    printf("tiempo de preparación de criba: %f segundos\n", segundosCriba);

	printf("Calculando Polinomio...\n");
	t_inicio = clock();
	crearMatrizNula(&qs_data);

	int res = 1;
	qs_data.intervalo.Qxi = NULL;

	long polinomio_count = 0;
	long prev_n_BSuaves = qs_data.n_BSuaves;
	long full_relations = 0;
	while(res==1){
		// MPQS: generar nuevo polinomio
		generate_mpqs_poly(&qs_data);

		// === CRIBA LOGARÍTMICA ===
		// Pre-filtrar candidatos con criba logarítmica antes de trial division
		long *sieve_candidates = NULL;
		unsigned long n_candidates = 0;
		sieve_mpqs(&qs_data, xmax, &sieve_candidates, &n_candidates);

		if (n_candidates == 0) {
			free(sieve_candidates);
			polinomio_count++;
			continue;
		}

		// Construir Xi y Qxi solo para los candidatos de la criba
		unsigned long npos = n_candidates;

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

		for (unsigned long i = 0; i < npos; i++) {
			mpz_init(qs_data.intervalo.Xi[i]);
			mpz_set_si(qs_data.intervalo.Xi[i], sieve_candidates[i]);
			mpz_init(qs_data.intervalo.Qxi[i]);
			eval_mpqs_Qx(&qs_data, qs_data.intervalo.Xi[i], qs_data.intervalo.Qxi[i]);
		}

		free(sieve_candidates);

		polinomio_count++;
		if (qs_data.blocks.length > 0) {
			res = factoringBlocks(&qs_data, npos, 0, xmax);
		} else {
			res = factoringTrial(&qs_data, npos, 0, xmax);
		}
		long found_this = qs_data.n_BSuaves - prev_n_BSuaves;
		if (found_this > 0) {
			/* Clasificar: full vs combined */
			full_relations += found_this; /* ajustado abajo */
			printf("Polinomio %ld: %lu candidatos criba → %ld B_suaves (total: %ld, parciales: %lu)\n",
				polinomio_count, npos, found_this, qs_data.n_BSuaves, qs_data.partials.n);
			fflush(stdout);
		}
		prev_n_BSuaves = qs_data.n_BSuaves;
	}
	printf("Polinomios procesados: %ld\n", polinomio_count);
	fflush(stdout);
	printf("Numeros B_Suaves encontrados:%ld\n",qs_data.n_BSuaves);
	printf("Parciales almacenadas sin emparejar: %lu\n", qs_data.partials.n);
	t_final = clock();
	double segundosPolinomio = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
	printf("tiempo de calculo del polinomio:%fs\n",segundosPolinomio);	
	
	
	printf("Escribiendo matriz...");
	imprimirMatriz(qs_data.mat);

	// Guardar roota en archivo para que mulPoli lo use
	FILE *fr = fopen("roota.txt", "w");
	if (fr) {
		mpz_out_str(fr, 10, qs_data.roota);
		fprintf(fr, "\n");
		fclose(fr);
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
	//TODO:cambiar bloques por punteros a base
	//reservo memoria para el array de bloques
	qs_data->blocks.block = (prime_block*)malloc((qs_data->blocks.length)*sizeof(prime_block));
	
	
	//reservo memoria para cada bloque
	for(int i = 0; i < qs_data->blocks.length; i++)
	{
		qs_data->blocks.block[i].factors = (prime*)malloc(n*sizeof(prime));//reservo memoria para n factores
		//printf("%x\n",qs_data->blocks.block[i].factors);
	}
		
	mpz_t mulTemp;//variable multiplicacion de bloques
	mpz_init(mulTemp);
	mpz_set_ui(mulTemp,1);
	
	//Creo bloques de tamaño n a partir de los factores de la base que
	//esta almacenada en la estructura
	int contBlock = 0;
	int contFact = 0;
	
	for (int i = 0; i < qs_data->base.length; i++)
	{
		//si el blo1ue se llena avanzo al siguiente
		if(contFact==n)
		{
			qs_data->blocks.block[contBlock].length = contFact;
			mpz_init(qs_data->blocks.block[contBlock].prod_factors);
			mpz_set(qs_data->blocks.block[contBlock].prod_factors,mulTemp);
			mpz_set_ui(mulTemp,1);
			
			contBlock++;
			contFact = 0;
		}
		//TODO:cambiar value de factors por puntero
		mpz_init(qs_data->blocks.block[contBlock].factors[contFact].value);
		mpz_set(qs_data->blocks.block[contBlock].factors[contFact].value,qs_data->base.primes[i].value);
		mpz_mul(mulTemp,mulTemp,qs_data->base.primes[i].value);
		//gmp_printf("%Zd,",qs_data->blocks.block[contBlock].factors[contFact].value);
		contFact++;
	}
	
	//asigno el tamaño del ultimo bloque y el la
	//multiplicacion de los facatores del ultimo bloque
	qs_data->blocks.block[contBlock].length = contFact;
	mpz_init(qs_data->blocks.block[contBlock].prod_factors);
	mpz_set(qs_data->blocks.block[contBlock].prod_factors,mulTemp);
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
    gmp_printf("Primo mas grande en la base: %Zd\n", primes[contRes-1].value);
    printf("Primos leidos del archivo: %ld, residuos cuadraticos: %ld de %ld requeridos\n",
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
	
	while ((c = getopt(argc, argv, "d:h:b:c:")) != -1){
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
	fprintf(stderr,"Uso: ./B_QSieve (-d | -h) <N> [-b <NBLOCKS>] [-c <NCORES>] \n");
	fprintf(stderr,"Opciones:\n");
	fprintf(stderr,"-d	# Especifica que el numero N es decimal\n");
	fprintf(stderr,"-h	# Especifica que el numero N es hexadecimal\n");
	fprintf(stderr,"-b	# Al usar esta opcion se deben especificar el numero de Bloques\n");
	fprintf(stderr,"-C	# Especifica el numero de procesadores logicos que se quieren usar\n");
}
