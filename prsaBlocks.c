/**
 * @file 
 * @author Jhon Alejandro Melo <alejandromelo@unicauca.edu.co>
 * 			Juan Manuel Campo <>
 * @copyright GNU Public License. 
 *
 * @brief Implementacion del algoritmo Quadratic Sieve
 */

#include <stdio.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "split.h"
#include "structsqs.h"
#include <unistd.h>

void agregarAVectorBlock();
void bloques(int n);
int buscarSuave(mpz_t suave);
int calcularResiduos(mpz_t n, mpz_t longitud);
void crearMatrizNula(int filas, int columnas);
void insertarNumero(int fila, int columna, int posFila);
void intervaloPolinomio(mpfr_t n, mpz_t result);
void imprimirMatriz(int filas, int columnas);
void longitudBase(mpfr_t n, mpz_t result);
void polinomioFermat(mpz_t n, mpz_t intervalo, int numSuaves);
int validarResultado(mpz_t result);

void reducirPolinomio();

//TODO:quitar variables globales
int **matriz;
int posFila = 0;//pos fila matriz para agregar 
int n_block;
int cPrimes;
//almacena informacion para recuperar el vector de los numeros suaves
data_div_table block_table;

int main(int argc, char * argv[]){	
	if(argc != 3){
		fprintf(stderr,"Parametros no validos\n");
		fprintf(stderr,"Debe especificar el numero que se va factorizar y la cantidad de primos por bloque\n");
		fprintf(stderr,"Ejemplo: ./prsa 15347 10\n");
		exit(EXIT_FAILURE);
	}
	
	n_block = atoi(argv[2]);//numeros por bloque 
	
	FILE * fgcTemp;
	fgcTemp = fopen("gcd.txt","w");
	fclose(fgcTemp);
	
	
	clock_t t_inicio, t_final;
	
	double segundos = 0;
	
	//declaracion de variables
	mpz_t n;
	mpz_t intervaloP;
	mpz_t longitudB;
	mpfr_t num;
	
	//instancia de variables
	mpz_init(n);
	mpz_init(intervaloP);
	mpz_init(longitudB);
	mpfr_init(num);
		
	//digitos de N
	printf("Numero de digitos: %lu\n",strlen(argv[1]));
	mpz_set_str(n, argv[1], 10);
	//mpz_out_str(stdout,10,n);
	//printf("\n");
	
	//longitud de la base de primos
	mpfr_set_str(num, argv[1], 10, MPFR_RNDU);
	longitudBase(num,longitudB);
	printf("Longitud de la base:");
	mpz_out_str(stdout,10,longitudB);
	printf("\n");
	
	cPrimes = mpz_get_ui(longitudB);//cantidad de primos
	if(cPrimes <= n_block)
		n_block = cPrimes;
	
	block_table.data = (data_div*)malloc((cPrimes+1)*sizeof(data_div));
	crearMatrizNula(cPrimes+1,cPrimes);
	
	
	//intervalo del polinomio
	mpfr_set_str(num, argv[1], 10, MPFR_RNDU);
	intervaloPolinomio(num,intervaloP);
	printf("intervalo del polinomio: ");
	mpz_out_str(stdout,10,intervaloP);
	printf("\n");
	
	t_inicio = clock();
	//generar residuos
	printf("Generando residuos cuadraticos...\n");
	int residuos = calcularResiduos(n,longitudB);
	printf("%d residuos generados en residuos.txt\n",residuos);
	t_final = clock();
	segundos = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
	printf("tiempo de calculo de residuos:%fs\n",segundos);

	printf("Creando bloques...\n");
	bloques(n_block);
	printf("Bloques creados: %.1f\n",ceil((float)cPrimes/n_block));
	fflush(stdout);

	//Calcular polinomio de Fermat por metodo de bloques
	t_inicio = clock();
	printf("Calculando polinomio de Fermat con reduccion por bloques\n");
	polinomioFermat(n,intervaloP,cPrimes+1);
	t_final = clock();
	segundos = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
	printf("tiempo de calculo del polinomio:%fs\n",segundos);
	
	imprimirMatriz(cPrimes+1,cPrimes);
	
	//Liberar memoria
	mpz_clear(n);
	mpz_clear(intervaloP);
	mpz_clear(longitudB);
	mpfr_clear(num);
	free(block_table.data);
	free(matriz);
	
	exit(EXIT_SUCCESS);
}


void crearMatrizNula(int filas, int columnas){
	//reservar memoria para matriz
	matriz = (int**)malloc(filas*sizeof(int*));
	
	for (int i = 0; i < filas; i++) 
	{
		matriz[i] = (int*)malloc(columnas*sizeof(int));
		memset(matriz[i],0,columnas*sizeof(int));
	}
}

void imprimirMatriz(int filas, int columnas){
	FILE* f = fopen("matrix.txt","w");
	printf("%d %d\n",filas, columnas);
	fprintf(f,"%d %d\n",filas,columnas);
	int v[columnas];
	int cont;
	for (int i = 0; i < filas; i++)
	{
		cont = 0;
		for (int j = 0; j < columnas; j++)
		{
			printf("%d ",matriz[i][j]);
			if(matriz[i][j]==1){
				v[cont] = j;
				cont++;
			}
		}
		fprintf(f,"%d ",cont);
		for (int k = 0; k < cont; k++)
		{
			fprintf(f,"%d ",v[k]);
		}
		printf("\n");	
		fprintf(f,"\n");
	}
	
	
	fclose(f);
}	

/**
 * @brief calcula los residuos cuadraticos del numero n usando la funcion legendre,
 * los primos se obtienen del archivo primes.txt y los almacena en un archivo residuos.txt 
 * @param n: numero que se le evaluaran los residuos cuadraticos
 * @param longitud: numero de residuos que se requieren para el numero n
 * @return retorna el numero de 
 * residuos encontrados
 */
int calcularResiduos(mpz_t n, mpz_t longitud){
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
			mpz_init_set_str(p, buf, 10);//se asigna un numero del archivo a p
			
			//si n es residuo cuadratico mod p se agrega al archivo
			if((mpz_legendre(n,p)==1)){
				contRes++;
				if((fr = fopen("residuos.txt","a"))==NULL){
					perror("fopen");
					exit(EXIT_FAILURE);
				}
				fprintf(fr,"%s",buf);
				fclose(fr);
			}else if(mpz_cmp_ui(p,2)==0){
				//si el primo es 2 es residuo cuadratico de todo n
				contRes++;
				if((fr = fopen("residuos.txt","a"))==NULL){
					perror("fopen");
					exit(EXIT_FAILURE);
				}
				fprintf(fr,"%s",buf);
				fclose(fr);
			}
			//si los residuos es igual a la longitud termina
			if((mpz_cmp_ui(longitud,contRes)==0)){
				mpz_clear(p);
				break;
			}
			mpz_clear(p);
		}
	}
	printf("Se usaron %d primos\n",contn);  
	return contRes;
}

/**
 * @brief calcula el numero de residuos que se necesitan
 * para factorizar el numero n
 * @param n: numero que se le calcula la longitud de la base
 * @param result: variable en la que se devuelve el valor calculado
 */
void longitudBase(mpfr_t n, mpz_t result){
	//formuala: result = ((e^sqrt(ln(n)*ln(ln(n))))^(sqrt(2)/4))
	mpfr_t ln1, ln2,e,pow;//variables
	
	//inicializacion de variables
	mpfr_init(ln1);
	mpfr_init(ln2);
	mpfr_init(e);
	mpfr_init(pow);
	
	mpfr_set_str(e, "2.71828182845904523536", 10, MPFR_RNDZ);//define euler
	mpfr_set_str(pow, "0.3535533905932738", 10, MPFR_RNDZ);//define sqrt(2)/4
	mpfr_log(ln1, n, MPFR_RNDZ);//ln1= log(n)
	mpfr_log(ln2, ln1,MPFR_RNDZ);//ln2=log(ln1)
	mpfr_mul(n, ln1, ln2, MPFR_RNDZ);//n=ln1*ln2
	mpfr_sqrt(n, n, MPFR_RNDZ);//sqrt(n)
	mpfr_pow (n, e, n, MPFR_RNDZ);//n=e^n
	mpfr_pow (n, n, pow, MPFR_RNDZ);//n=n^pow
	mpfr_get_z(result, n, MPFR_RNDZ);//se le asigna a result la parte entera de n
	
	mpfr_clear(ln1);
	mpfr_clear(ln2);
	mpfr_clear(e);
	mpfr_clear(pow);
}

/**
 * @brief calcula el intervalo en el que se deben probar los residuos 
 * cuadraticos del numero n
 * @param n: numero que se le calcula el intervalo que necesita
 * @param result: variable que se ultiliza para devolver el resultado
 */
void intervaloPolinomio(mpfr_t n, mpz_t result){
	//formuala: result = ((e^sqrt(ln(n)*ln(ln(n))))^(sqrt(2)/4))^3
	mpfr_t ln1, ln2,e,pow;
	
	mpfr_init(ln1);
	mpfr_init(ln2);
	mpfr_init(e);
	mpfr_init(pow);
	
	mpfr_set_str(e, "2.71828182845904523536", 10, MPFR_RNDZ);//define euler
	mpfr_set_str(pow, "0.3535533905932738", 10, MPFR_RNDZ);//define sqrt(2)/4
	mpfr_log(ln1, n, MPFR_RNDZ);//ln1 = log(n)
	mpfr_log(ln2, ln1,MPFR_RNDZ);//ln2 = log(log(n))
	mpfr_mul(n, ln1, ln2, MPFR_RNDZ);//n = ln1*ln1
	mpfr_sqrt(n, n, MPFR_RNDZ);// n = sqrt(n)
	mpfr_pow (n, e, n, MPFR_RNDZ);// n = e^n
	mpfr_pow (n, n, pow, MPFR_RNDZ);// n = n^pow
	mpfr_pow_si(n, n, 3, MPFR_RNDZ);//n = n^3
	//mpfr_out_str(stdout, 10, 0, n, MPFR_RNDZ);
	mpfr_get_z(result, n, MPFR_RNDZ);//se le asigna a result la parte entera de n
	mpfr_clear(ln1);
	mpfr_clear(ln2);
	mpfr_clear(e);
	mpfr_clear(pow);
}

/**
 * @brief Calcula los valores del polinomio y los agrega al archivo polinomio.txt con el metodo de los bloques
 * @param n: Numero al que se le calculan los valores del polinomio
 * @param intervalo: Nos indica el rango en el que deben evaluarse los residuos del numero n
 */
void polinomioFermat(mpz_t n, mpz_t intervalo, int numSuaves){
	//Q(Xi)=(sqrt(n)+i)^2-n
	mpz_t raiz;//sqrt(n)
	mpz_t limite;//limite del ciclo
	mpz_t result;//resultados del polinomio
	
	mpz_init(raiz);
	mpz_init(limite);
	mpz_init(result);
	
	mpz_set(limite,intervalo);//limite = intervalo
	mpz_set_ui(intervalo,0);
	mpz_sqrt(raiz,n);//raiz = sqrt(n)
	//mpz_mul_si(intervalo,intervalo,-1);//intervalo = -1*intervalo
	
	FILE * fp;
	if((fp = fopen("polinomioB.txt","w"))==NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	fclose(fp);

	char buf[BUFSIZ];
	
	while(mpz_cmp(intervalo,limite)!=0){
		mpz_mul_si(intervalo,intervalo,-1);
		if(numSuaves == 0)break;//condicion para saber si se encontraron todos los suaves que se necesitan
		//se incrementa el valor del intervalo en 1
		if(mpz_sgn(intervalo)==1 || mpz_sgn(intervalo)==0)
			mpz_add_ui(intervalo,intervalo,1);
		
		//calcular resultado
		mpz_add(result,raiz,intervalo);//result=raiz+intervalo
		mpz_pow_ui(result,result,2);//result = result^2
		mpz_sub(result,result,n);//result = result-n
		
		//Si al validar el resultado da 1 se agrega al archivo polinomio
		if(validarResultado(result)){
			if(buscarSuave(result))continue;//si el numero del polinomio ya esta continua con el siguiente
			numSuaves--;//se encontro un suave no repetido
			agregarAVectorBlock();		
			posFila++;
			//agregar al archivo
			memset(buf,0,BUFSIZ);
			if((fp = fopen("polinomioB.txt","a"))==NULL){
				perror("fopen");
				exit(EXIT_FAILURE);
			}
			mpz_t Xi;
			mpz_init(Xi);
			mpz_add(Xi,raiz,intervalo);
			mpz_out_str(fp,10,Xi);
			fprintf(fp,";");
			mpz_out_str(fp,10,result);
			fprintf(fp,"\n");
			fclose(fp);
		}
	}
	
	mpz_clear(raiz);
	mpz_clear(limite);
	mpz_clear(result);
}

/**
 * @brief Busca el numero suave si ya esta almacenado en el archivo
 * @param suave: numero suave que sale del polinomio
 * @return retorna 0 si el numero no esta en el archivo o 1 si el numero ya esta almacenado
 */
int buscarSuave(mpz_t suave){
	char buf[BUFSIZ];
	FILE * fp;
	//validar si el numero ya esta
	if((fp = fopen("polinomioB.txt","r"))==NULL){
			perror("fopen");
			exit(EXIT_FAILURE);
	}
	mpz_t numTemp;
	int bd = 0;
	while(!feof(fp)){
		memset(buf,0,BUFSIZ);
		if(fgets(buf,BUFSIZ,fp)!=NULL){
			
			char ** tokens;
			int n = 0;
			tokens = split(buf,";",&n);
			if(n==0)break;
			mpz_init(numTemp);
			mpz_set_str(numTemp,tokens[1],10);
			free(tokens);
			if(mpz_cmp(suave,numTemp)==0){
				//printf("Colision");
				bd=1;
				mpz_clear(numTemp);
				break;
			}
			mpz_clear(numTemp);
		}
	}
	fclose(fp);
	
	return bd;
}

/**
 * @brief crea un archivo bloques.txt donde se almacena la multiplicacion
 * de los residuos separados en bloques de tamaño n
 * @param n:tamaño de los bloques
 */
void bloques(int n){
	FILE * fb;
	fb = fopen("bloques.txt","w");
	
	FILE * fbf;
	fbf = fopen("bloquesFac.txt","w");
	
	fclose(fbf);
	fclose(fb);

	char buf[BUFSIZ];

	int tb=n;//tamaño de bloques

	mpz_t mulTemp;//variable multiplicacion de bloques
	mpz_init(mulTemp);
	mpz_set_ui(mulTemp,1);

	mpz_t valorTemp;//variable para almacenar el numero del archivo

	FILE * fr;
	fr = fopen("residuos.txt","r");

	while(!feof(fr)){
		while(tb!=0 && !feof(fr)){
			memset(buf,0,BUFSIZ);
			if(fgets(buf,BUFSIZ,fr)!=NULL){
				buf[strlen(buf)-1] = '\0';
				fbf = fopen("bloquesFac.txt","a");
				fprintf(fbf,"%s;",buf);
				fclose(fbf);
				mpz_init_set_str(valorTemp, buf, 10);
				mpz_mul(mulTemp,mulTemp,valorTemp);//mulTemp=mulTemp*valorTemp
				tb--;
				mpz_clear(valorTemp);
  			}
  		}
  		fbf = fopen("bloquesFac.txt","a");
		fprintf(fbf,"\n");
		fclose(fbf);
  		tb=n;
  		memset(buf,0,BUFSIZ);
		fb = fopen("bloques.txt","a");
		mpz_out_str(fb,10,mulTemp);
		fprintf(fb,"\n");
		fclose(fb);
		mpz_set_ui(mulTemp,1);
	}
	fclose(fr);
	mpz_clear(mulTemp);
}

/**
 * @brief Valida si un numero del polinomio se divide en la base de residuos
 * usando el metodo de los bloques, si el numero es liso se almacenan sus factores en un 
 * archivo
 * @param result: Numero que se valida si es divisible en la base
 * @return retorna 1 si es divisible en la base o 0 si no lo es
 */
 //TODO:Cambiar nombre parametro
int validarResultado(mpz_t result){
	
	char buf[BUFSIZ];
	mpz_t resultf;
	mpz_init(resultf);//variable de resultado final
	mpz_set(resultf,result);//resultf = result
	//si el resultado es negativo se vuelve positivo
	if(mpz_sgn(resultf)==-1)
		mpz_mul_si(resultf,resultf,-1);

	mpz_t valorB;//valor del bloque

	mpz_t mcd;//valor del mcd
	mpz_init(mcd);
	
	memset(block_table.data,0,(cPrimes+1)*sizeof(data_div));
	FILE * fb;
	fb = fopen("bloques.txt","r");
	int contBloque = 0;
	int i = 0;
	while(!feof(fb)){
		
		memset(buf,0,BUFSIZ);
		
		if(fgets(buf,BUFSIZ,fb)==NULL){
			break;
		}
		
		contBloque++;
		
		mpz_init_set_str(valorB, buf, 10);
		mpz_t mcdA;//variable para controlar cuando cambia el valor del mcd
		mpz_init(mcdA);
		
		//se le asigna el valor de mcd a mcdA
		mpz_gcd(mcd,resultf,valorB);
		mpz_set(mcdA,mcd);
		int contGcda = 0;
		
		//TODO:Optimizar ciclo
		do{
			mpz_gcd(mcd,resultf,valorB);
			//si el maximo comun divisor cambia se guarda en un archivo temporal
			if(mpz_cmp(mcd,mcdA)!=0){
				//Se almacena en una estructura el gcd
				//las veces que se repite y el bloque al que pertenece
				block_table.data[i].block = contBloque;
				block_table.data[i].gcd = mpz_get_ui(mcdA);
				block_table.data[i].periodo = contGcda;
				mpz_set(mcdA,mcd);
				contGcda = 1;
				i++;
			}else{
				contGcda++;
			}
			//si el maximo comun divisor es 1 sale del ciclo
			if(mpz_cmp_ui(mcd,1)==0){
				break;
			}
			
			mpz_divexact(resultf,resultf,mcd);
		}while(mpz_cmp_ui(mcd,1)!=0);
		
		mpz_clear(valorB);
		mpz_clear(mcdA);
	}
	fclose(fb);
	mpz_clear(mcd);
	if(mpz_cmp_ui(resultf,1)==0){
		mpz_clear(resultf);
		//mpz_out_str(stdout,10,result);
		//printf("\n");
		return 1;
	}
	mpz_clear(resultf);
	return 0;
}

/**
 * @brief 
 * @param 
 * @return
 */
void agregarAVectorBlock(){
	//TODO:Cambiar estrucura y quitar lectura del archivo, la estructura ya guarda bien la informacion
	int i = 0;
	/*if(posFila==105){
		while(block_table.data[i].gcd!=0){
			int contV = block_table.data[i].periodo;//cliclos que se repite el mcd
			int block = block_table.data[i].block;//bloque en el que esta el mcd
			int mcdE = block_table.data[i].gcd;
			i++;
			printf("periodo:%d, bloque:%d, gcd:%d\n",contV,block,mcdE);
		}
	}
	i=0;*/
	while(block_table.data[i].gcd!=0){
		/*if(posFila==105){
			printf("i=%d\n",i);	
			printf("%d,%d,%d\n",block_table.data[i].gcd,block_table.data[i].periodo,block_table.data[i].block);
		}*/
		FILE * fblocks;
		fblocks = fopen("bloquesFac.txt","r");
	
		mpz_t mcd;//valor del mcd
		mpz_init(mcd);
		mpz_set_ui(mcd,block_table.data[i].gcd);
		
		int contV = block_table.data[i].periodo;//cliclos que se repite el mcd
		int block = block_table.data[i].block;//bloque en el que esta el mcd
		//int mcdE = block_table.data[i].gcd;
		i++;
		//TODO: arreglar ciclo
		//si los periodos son pares se ingresan ceros
		if((contV%2)==0){
			//printf("fila:%d Block:%d mcd:%d par\n",posFila,block,mcdE);
			for (int i = 0; i < n_block; i++)
			{
				insertarNumero(posFila,(block-1)*n_block+i,0);
			}
			continue;
		}
		char buf2[BUFSIZ];
		int contBlock = 0;
		while(!feof(fblocks)){
			memset(buf2,0,BUFSIZ);
			if(fgets(buf2,BUFSIZ,fblocks)==NULL){
				break;
			}
			contBlock++;
			if(contBlock == block){
				char ** tokens2;
				int n2 = 0;
				tokens2 = split(buf2,";",&n2);
				
				for (int i = 0; i < n2-1; i++)
				{
					
					int nblock = atoi(tokens2[i]);
					/*if(posFila==105){
						mpz_out_str(stdout,10,mcd);
						printf(";P=%d\n",nblock);
						fflush(stdout);
					}*/
					if(mpz_divisible_ui_p(mcd,nblock)){
						insertarNumero(posFila,(block-1)*n_block+i,contV%2);
					}else{
						insertarNumero(posFila,(block-1)*n_block+i,0);
					}
				}
				free(tokens2);
			}
		}
		
		//printf("%d; %d\n",contV,block);
		
		mpz_clear(mcd);
		fclose(fblocks);
		/*if(posFila==105){
			for (int i = 0; i < 115; i++)
			{
				printf("%d ",matriz[105][i]);
			}
		}*/
	}
}

void insertarNumero(int fila, int columna, int valor){
	
	/*if(fila==105){
		printf("Fila:%d Columna:%d Valor:%d\n",fila,columna,valor);
		printf("Fila:%d Columna:%d Valor:%d\n",fila,columna,matriz[fila][columna]);
	}
	fflush(stdout);*/
	if(matriz[fila][columna] == 0 && valor == 0){
		matriz[fila][columna] = 0;
		return;
	}
	
	if(matriz[fila][columna] == 1 && valor == 1){
		matriz[fila][columna] = 0;
		return;
	}
	
	if(matriz[fila][columna] == 0 && valor == 1){
		matriz[fila][columna] = 1;
		return;
	}
	
	if(matriz[fila][columna] == 1 && valor == 0){
		matriz[fila][columna] = 1;
		return;
	}
}
