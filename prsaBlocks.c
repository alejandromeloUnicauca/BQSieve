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
int calcularResiduos(mpz_t n, mpz_t longitud);
void crearMatrizNula(int filas, int columnas);
void insertarNumero(int fila, int columna, int posFila);
void intervaloPolinomio(mpfr_t n, mpz_t result);
void imprimirMatriz(int filas, int columnas);
void longitudBase(mpfr_t n, mpz_t result);
void polinomioFermat(mpz_t n, mpz_t intervalo, int numSuaves);
int validarResultado(mpz_t result);

void reducirPolinomio();


int **matriz;
int posFila = 0;//pos fila matriz para agregar 
int n_block;
//almacena informacion para recuperar el vector de los numeros suaves
blocks_div_table * block_table;
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
	mpz_out_str(stdout,10,n);
	printf("\n");
	
	//longitud de la base de primos
	mpfr_set_str(num, argv[1], 10, MPFR_RNDU);
	longitudBase(num,longitudB);
	printf("Longitud de la base:");
	mpz_out_str(stdout,10,longitudB);
	printf("\n");
	
	int cPrimes = mpz_get_ui(longitudB);//cantidad de primos
	if(cPrimes <= n_block)
		n_block = cPrimes;
	
	
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
	
	for (int i = 0; i < columnas; i++)
	{
		for (int j = 0; j < filas; j++)
		{
			printf("%d ",matriz[j][i]);
		}
		printf("\n");
		
	}
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
			if((mpz_cmp_ui(longitud,contRes)==0))
				break;
		}
	}
	
	mpz_clear(p);
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
	mpz_init(result);
	
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
			//validar si el numero ya esta
			if((fp = fopen("polinomioB.txt","r"))==NULL){
					perror("fopen");
					exit(EXIT_FAILURE);
			}
			int bd = 0;
			while(!feof(fp)){
				memset(buf,0,BUFSIZ);
				if(fgets(buf,BUFSIZ,fp)!=NULL){
					mpz_t numTemp;
					mpz_init(numTemp);
					mpz_set_str(numTemp,buf,10);
					if(mpz_cmp(result,numTemp)==0){
						//printf("Colision");
						bd=1;
						break;
					}
				}
			}
			fclose(fp);
			if(bd==1)continue;//si el numero del polinomio ya esta continua con el siguiente
			numSuaves--;//se encontro un suave no repetido
			agregarAVectorBlock();		
			posFila++;
			//agregar al archivo
			memset(buf,0,BUFSIZ);
			if((fp = fopen("polinomioB.txt","a"))==NULL){
				perror("fopen");
				exit(EXIT_FAILURE);
			}
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
	mpz_init(valorTemp);

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
	mpz_clear(valorTemp);
}

/**
 * @brief Valida si un numero del polinomio se divide en la base de residuos
 * usando el metodo de los bloques, si el numero es liso se almacenan sus factores en un 
 * archivo
 * @param result: Numero que se valida si es divisible en la base
 * @return retorna 1 si es divisible en la base o 0 si no lo es
 */
int validarResultado(mpz_t result){

	char buf[BUFSIZ];
	mpz_t resultf;
	mpz_init(resultf);//variable de resultado final
	mpz_set(resultf,result);//resultf = result
	//si el resultado es negativo se vuelve positivo
	if(mpz_sgn(resultf)==-1)
			mpz_mul_si(resultf,resultf,-1);

	mpz_t valorB;//valor del bloque
	mpz_init(valorB);

	mpz_t mcd;//valor del mcd
	mpz_init(mcd);
	
	
	//block_table = (blocks_div_table*)malloc(cPrimes*sizeof(blocks_div_table*));
	//TODO:Cambiar archivo temporal por demora
	FILE * fgcTemp;
	fgcTemp = fopen("gcdTemp.txt","w");
	fclose(fgcTemp);
	FILE * fb;
	fb = fopen("bloques.txt","r");
	int contBloque = 0;
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
		
		do{
			mpz_gcd(mcd,resultf,valorB);
			//si el maximo comun divisor cambia se guarda en un archivo temporal
			if(mpz_cmp(mcd,mcdA)==0){
				contGcda++;
			}else{
				//Se almacena en un archivo temporal el gcd
				//las veces que se repite y el bloque al que pertenece
				fgcTemp = fopen("gcdTemp.txt","a");
				mpz_out_str(fgcTemp,10,mcdA);		
				fprintf(fgcTemp,";%d",contGcda);
				fprintf(fgcTemp,";%d\n",contBloque);
				fclose(fgcTemp);
				mpz_set(mcdA,mcd);
				contGcda = 1;
			}
			//si el maximo comun divisor es 1 sale del ciclo
			if(mpz_cmp_ui(mcd,1)==0){
				break;
			}
			
			mpz_divexact(resultf,resultf,mcd);
		}while(mpz_cmp_ui(mcd,1)!=0);
		mpz_clear(mcdA);
	}
	fclose(fb);
	
	if(mpz_cmp_ui(resultf,1)==0){
		FILE * fgc = fopen("gcd.txt","a");
		fgcTemp = fopen("gcdTemp.txt","r");
		while(!feof(fgcTemp)){
			memset(buf,0,BUFSIZ);
			if(fgets(buf,BUFSIZ,fgcTemp)==NULL){
				break;
			}	
			fprintf(fgc,"%s",buf);
		}
		fclose(fgcTemp);
		fclose(fgc);
		
		return 1;
	}

	return 0;
}

void agregarAVectorBlock(){
	char buf[BUFSIZ];
	FILE * fgcdTemp;
	fgcdTemp = fopen("gcdTemp.txt","r");
	
	while(!feof(fgcdTemp)){
		memset(buf,0,BUFSIZ);
		if(fgets(buf,BUFSIZ,fgcdTemp)==NULL){
			break;
		}
		
		FILE * fblocks;
		fblocks = fopen("bloquesFac.txt","r");
	
		char ** tokens;
		int n = 0;
		tokens = split(buf,";",&n);
		mpz_t mcd;//valor del mcd
		mpz_init(mcd);
		mpz_set_str(mcd,tokens[0],10);
		
		int contV = atoi(tokens[1]);//cliclos que se repite el mcd
		int block = atoi(tokens[2]);//bloque en el que esta el mcd
		int mcdE = atoi(tokens[0]);
		//TODO: arreglar ciclo
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
					if(mpz_divisible_ui_p(mcd,nblock)){
						insertarNumero(posFila,(block-1)*n_block+i,contV%2);
					}else{
						insertarNumero(posFila,(block-1)*n_block+i,0);
					}
				}
			}
		}
		
		//printf("%d; %d\n",contV,block);
		
		/*TODO:leer archivo de bloquesFac
		 * dividir el mcd entre los numeros del bloque
		 * agregar los resultados al vector
		*/
		mpz_clear(mcd);
		fclose(fblocks);
	}
	fclose(fgcdTemp);	
}

void insertarNumero(int fila, int columna, int valor){
	
	//printf("Fila:%d Columna:%d Valor:%d\n",fila,columna,valor);
	
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
