#include <ctype.h>
#include <gmp.h>
#include <getopt.h>
#include <math.h>
#include <mpfr.h>
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
void getPrimesBaseLength(mpz_t n, long * result);
void getIntervalLength(mpz_t n, mpz_t result);
int generatePrimesBase(mpz_t n, long base_length, prime * primes);
void freeStruct(qs_struct * qs_data);
void usage();

int main(int argc, char **argv)
{
	int flagd = 0; 
	int flagh = 0;
	char *hdvalue = NULL;
	char *bvalue = NULL;
	int c;
	
	opterr = 0;
	
	while ((c = getopt(argc, argv, "d:h:b:")) != -1){
		switch(c){
			case 'd':
				if(flagh == 1){
					fprintf (stderr, "Solo puedes usar -h o -d pero no ambos\n");
					usage();
					exit(EXIT_FAILURE);
				}
				flagd = 1;
				hdvalue = optarg;
				break;
			case 'h':
				if(flagd == 1){
					fprintf (stderr, "Solo puedes usar -h o -d pero no ambos\n");
					usage();
					exit(EXIT_FAILURE);
				}
				flagh = 1;
				hdvalue = optarg;
				break;
			case 'b':
				bvalue = optarg;
				break;
			case '?':
				if (strchr("h", optopt) != NULL)
					fprintf (stderr, "La opción -%c requiere un argumento.\n", optopt);
				else if (strchr("d", optopt) != NULL)
					fprintf (stderr, "La opción -%c requiere un argumento.\n", optopt);
				else if (strchr("b", optopt) != NULL)
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
	
	//Declaracion de variables
	qs_struct qs_data;
	clock_t t_inicio, t_final;
	
	//Instancia de variables
	qs_data.n_BSuaves = 0;   
	mpz_inits(qs_data.n,qs_data.intervalo.length,NULL);
	if(bvalue!=NULL)qs_data.blocks.length = atol(bvalue);
	else qs_data.blocks.length = 0;
	
	
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
	
	//longitud de la base de primos
	getPrimesBaseLength(qs_data.n,&qs_data.base.length);
	printf("Longitud de la base de primos:%ld\n", qs_data.base.length);
	
	//Intervalo del polinomio
	getIntervalLength(qs_data.n,qs_data.intervalo.length);
	gmp_printf("Intervalo del polinomio:%Zd\n", qs_data.intervalo.length);
	
	//Generar base de primos
	qs_data.base.primes = (prime*)malloc((qs_data.base.length)*sizeof(prime));
	
	t_inicio = clock();
	printf("Generando base de primos...\n");
	int residuos = generatePrimesBase(qs_data.n,qs_data.base.length,qs_data.base.primes);
	t_final = clock();
	printf("Base de primos generada. %d primos en la base\n",residuos);
	
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
	
	long lengthXi = 0;
	
	printf("Cribando...\n");
	t_inicio = clock();
	long *Xi = sieving(&qs_data,&lengthXi); 
	t_final = clock();
	double segundosCriba = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
	printf("tiempo de cribado:%fs\n",segundosCriba);	
	
	
	qs_data.intervalo.length_Xi = lengthXi;
	qs_data.intervalo.length_Qxi = qs_data.intervalo.length_Xi;
	qs_data.intervalo.Xi = Xi;
	
	printf("Calculando Polinomio...\n");
	t_inicio = clock();
	unsigned long endPos = 0;
	long numLote = 1;
	//fermat(&qs_data);
	crearMatrizNula(&qs_data);

	int res = 1;
	unsigned long posXi = 0;
	unsigned long sizeLote = 150000;
	while(endPos<lengthXi && res==1){
		endPos=fermat2(&qs_data,numLote++,sizeLote);
		if(qs_data.blocks.length > 0){
			res = factoringBlocks(&qs_data,sizeLote,posXi*sizeLote);
		}
		else{
			res = factoringTrial(&qs_data,sizeLote,posXi*sizeLote);  
		}
		posXi++;
	}
	 
	
	t_final = clock();
	double segundosPolinomio = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
	printf("tiempo de calculo del polinomio:%fs\n",segundosPolinomio);	
	
	
	
	/*printf("Factorizando...\n");
	t_inicio = clock();
	if(qs_data.blocks.length > 0){
		factoringBlocks(&qs_data);
	}
	else{
		factoringTrial(&qs_data,0,qs_data.intervalo.length_Qxi);  
	}
	t_final = clock();
	double segundosFactorizacion = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
	printf("tiempo de factorizacion:%fs\n",segundosFactorizacion);	
	
	printf("tiempo de Total:%fs\n",segundosCriba+segundosPolinomio+segundosFactorizacion);	*/
	
	printf("Escribiendo matriz...");
	imprimirMatriz(qs_data.mat);
	
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
	qs_data->mat.n_rows = qs_data->base.length+1;
	qs_data->mat.n_cols = qs_data->base.length;
	
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
	//multiplicacion de los facotores del ultimo bloque
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
	//TODO:validar si se crearon bloques
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
	
	//liberar memoria del intervalo
	/*mpz_clears(qs_data->n,qs_data->intervalo.length,NULL);
	
	for (long i = 0; i < qs_data->intervalo.length_Qxi; i++)
	{
		mpz_clear(qs_data->intervalo.Qxi[i]);
	}*/
	 
	free(qs_data->intervalo.Qxi);
	free(qs_data->intervalo.Xi);
	
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
int generatePrimesBase(mpz_t n, long base_length, prime * primes){
	//TODO:Quitar archivo de residuos.txt
	int contRes = 0;//contador de residuos encontrados

	mpz_t p;//variable temporal para los primos del archivo
	mpz_init(p);
	
	FILE * file;//file primes
	FILE * fr;//file residuos
	if((fr = fopen("residuos.txt","w")) == NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	fclose(fr);
	
	int contn = 0;//contador de primos que se ultilizan de primes.txt
	
	//si el archivo primes.txt no existe termina
	if ((file = fopen("primes.txt", "r")) == NULL) // open file
	{
		fprintf(stderr,"Falta archivo primes.txt");
		exit(EXIT_FAILURE);
	}
	
	char buf[BUFSIZ];
	while(!feof(file)){
		contn++;
		memset(buf, 0, BUFSIZ);
		if(fgets(buf,BUFSIZ,file)!=NULL){
			mpz_set_str(p, buf, 10);//se asigna un numero del archivo a p
			
			//si n es residuo cuadratico mod p se agrega al archivo
			if((mpz_legendre(n,p)==1) || (mpz_cmp_ui(p,2)==0)){
				//asigno memoria a los valores de prime
				mpz_init(primes[contRes].value);
				mpfr_init(primes[contRes].log_value);
				
				//almaceno el primo y el logaritmo del primo
				mpz_set(primes[contRes].value,p);
				
				
				mpfr_t pTemp;
				mpfr_init(pTemp);
				mpfr_set_z(pTemp,p,MPFR_RNDZ);
				mpfr_log(primes[contRes].log_value, pTemp, MPFR_RNDZ);//ln(p)
				primes[contRes].llog_value = mpfr_get_ui(primes[contRes].log_value,MPFR_RNDZ);
				
				mpfr_clear(pTemp);
				if((fr = fopen("residuos.txt","a"))==NULL){
					perror("fopen");
					exit(EXIT_FAILURE);
				}
				
				contRes++;
				
				fprintf(fr,"%s",buf);
				fclose(fr);
			}
			
			//si los residuos es igual a la longitud termina
			if((base_length==contRes)){
				break;
			}
			
		}
	}
	mpz_clear(p);
	fclose(file);
	printf("Se usaron %d primos\n",contn);  
	return contRes;
}

/**
 * @brief calcula el numero de residuos que se necesitan
 * para factorizar el numero n
 * @param n: numero que se le calcula la longitud de la base
 * @param result: variable en la que se devuelve el valor calculado
 */
void getPrimesBaseLength(mpz_t n, long * result){
	//formuala: result = ((e^sqrt(ln(n)*ln(ln(n))))^(sqrt(2)/4))
	
	//Declaracion de variables
	mpfr_t num, ln1, ln2, e, pow;
	
	//inicializacion de variables
	mpfr_inits(num,ln1,ln2,e,pow,NULL);

	mpfr_set_z(num,n,MPFR_RNDZ);
	mpfr_set_str(e, "2.71828182845904523536", 10, MPFR_RNDZ);//define euler
	mpfr_set_str(pow, "0.3535533905932738", 10, MPFR_RNDZ);//define sqrt(2)/4
	mpfr_log(ln1, num, MPFR_RNDZ);//ln1= log(num)
	mpfr_log(ln2, ln1,MPFR_RNDZ);//ln2=log(ln1)
	mpfr_mul(num, ln1, ln2, MPFR_RNDZ);//num=ln1*ln2
	mpfr_sqrt(num, num, MPFR_RNDZ);//sqrt(num)
	mpfr_pow (num, e, num, MPFR_RNDZ);//num=e^n
	mpfr_pow (num, num, pow, MPFR_RNDZ);//num=num^pow
	//mpfr_get_z(result, num, MPFR_RNDZ);//se le asigna a result la parte entera de num
	
	*result = mpfr_get_ui(num,MPFR_RNDZ);
	mpfr_clears(num,ln1,ln2,e,pow,NULL);
}

/**
 * @brief calcula el intervalo en el que se deben probar los residuos 
 * cuadraticos del numero n
 * @param n: numero que se le calcula el intervalo que necesita
 * @param result: variable que se ultiliza para devolver el resultado
 */
void getIntervalLength(mpz_t n, mpz_t result){
	//formuala: result = ((e^sqrt(ln(n)*ln(ln(n))))^(sqrt(2)/4))^3
	
	//Declaracion de variables
	mpfr_t num, ln1, ln2,e,pow;
	
	//inicializacion de variables
	mpfr_inits(num,ln1,ln2,e,pow,NULL);
	
	mpfr_set_z(num,n,MPFR_RNDZ);
	mpfr_set_str(e, "2.71828182845904523536", 10, MPFR_RNDZ);//define euler
	mpfr_set_str(pow, "0.3535533905932738", 10, MPFR_RNDZ);//define sqrt(2)/4
	mpfr_log(ln1, num, MPFR_RNDZ);//ln1 = log(num)
	mpfr_log(ln2, ln1,MPFR_RNDZ);//ln2 = log(log(num))
	mpfr_mul(num, ln1, ln2, MPFR_RNDZ);//num = ln1*ln1
	mpfr_sqrt(num, num, MPFR_RNDZ);// num = sqrt(num)
	mpfr_pow (num, e, num, MPFR_RNDZ);// num = e^num
	mpfr_pow (num, num, pow, MPFR_RNDZ);// num = num^pow
	mpfr_pow_si(num, num, 3, MPFR_RNDZ);//num = num^3
	mpfr_get_z(result, num, MPFR_RNDZ);//se le asigna a result la parte entera de num
	
	//Liberar memoria
	mpfr_clears(ln1,ln2,e,pow,num,NULL);
}

void usage(){
	fprintf(stderr,"Uso: ./B_QSieve (-d | -h) <N> [-b <NBLOCKS>] \n");
	fprintf(stderr,"Opciones:\n");
	fprintf(stderr,"-d	# Especifica que el numero N es decimal\n");
	fprintf(stderr,"-h	# Especifica que el numero N es hexadecimal\n");
	fprintf(stderr,"-b	# Al usar esta opcion se deben especificar el numero de Bloques\n");
}
