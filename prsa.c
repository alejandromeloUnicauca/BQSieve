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
#include <stdlib.h>
#include <string.h>
#include "split.h"
#include "structsqs.h"
#include <unistd.h>

void agregarAVectorDiv();
int calcularResiduos(mpz_t n, mpz_t longitud);
void crearMatrizNula(int filas, int columnas);
void insertarNumero(int fila, int columna, int posFila);
void intervaloPolinomio(mpfr_t n, mpz_t result);
void imprimirMatriz(int filas, int columnas);
void longitudBase(mpfr_t n, mpz_t result);
void polinomioFermat(mpz_t n, mpz_t intervalo, int numSuaves);
int validarResultado(mpz_t result);


int **matriz;
int posFila = 0;//pos fila matriz
int i_long;
data_divT *data_d;

int main(int argc, char * argv[]){	
	//TODO: obtimizar y comentar
	if(argc != 2){
		fprintf(stderr,"Parametros no validos\n");
		fprintf(stderr,"Debe especificar el numero que se va factorizar");
		exit(EXIT_FAILURE);
	}
	
	clock_t t_inicio, t_final;
	
	double segundos = 0;
	
	//declaracion de variables
	mpz_t n;
	mpz_t intervaloP;
	mpz_t longitudB;
	mpfr_t num;
	
	//instancia de variables
	mpz_init(n);
	mpz_init(longitudB);
	mpfr_init(num);
		
	//digitos de N
	printf("Numero de digitos: %lu\n",strlen(argv[1]));
	mpz_set_str(n, argv[1], 10);
	
	//longitud de la base de primos
	mpfr_set_str(num, argv[1], 10, MPFR_RNDU);
	longitudBase(num,longitudB);
	printf("Longitud de la base:");
	mpz_out_str(stdout,10,longitudB);
	printf("\n");
	
	i_long = mpz_get_ui(longitudB);//cantidad de primos
	data_d = (data_divT*)malloc((i_long+1)*sizeof(data_divT));
	
	crearMatrizNula(i_long+1,i_long);
	//imprimirMatriz(i_long+10,i_long);
	
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

	
	//crear los bloques
	
	//Calcular polinomio de Fermat
	t_inicio = clock();
	printf("Calculando polinomio de Fermat\n");
	polinomioFermat(n,intervaloP,i_long+1);
	t_final = clock();
	segundos = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
	printf("tiempo de calculo del polinomio:%fs\n",segundos);
	
	
	imprimirMatriz(i_long+1,i_long);
	
	//Liberar espacio de variables
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
			mpz_clear(p);
			//si los residuos es igual a la longitud termina
			if((mpz_cmp_ui(longitud,contRes)==0))
				break;
		}
	}
	printf("Se usaron %d primos\n",contn);  
	return contRes;
}

/***
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
 * @brief Calcula los valores del polinomio y los agrega al archivo polinomio.txt
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
	if((fp = fopen("polinomio.txt","w"))==NULL){
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
			if((fp = fopen("polinomio.txt","r"))==NULL){
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
					if(mpz_cmp(result,numTemp)==0){printf("Colision");bd=1;break;}
					mpz_clear(numTemp);
				}
			}
			fclose(fp);
			
			if(bd==1)continue;//si el numero del polinomio ya esta continua con el siguiente
			numSuaves--;//se encontro un suave no repetido
			agregarAVectorDiv();		
			posFila++;
			//agregar al archivo
			memset(buf,0,BUFSIZ);
			if((fp = fopen("polinomio.txt","a"))==NULL){
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
			mpz_clear(Xi);
		}
	}
	
	mpz_clear(raiz);
	mpz_clear(limite);
	mpz_clear(result);
}

/**
 * @brief Valida si un numero del polinomio se divide en la base de residuos
 * @param result: Numero que se valida si es divisible en la base
 * @return retorna 1 si es divisible en la base o 0 si no lo es
 */
int validarResultado(mpz_t result){
	//mpz_out_str(stdout,10,result);
	//printf("\n");
	mpz_t p;//primo
	mpz_t resultf;
	
	mpz_init(resultf);//variable de resultado final
	mpz_set(resultf,result);//resultf = result
	//si el resultado es negativo se vuelve positivo
	if(mpz_sgn(resultf)==-1)
		mpz_mul_si(resultf,resultf,-1);
		
	memset(data_d,0,(i_long+1)*sizeof(data_divT));
	FILE * fr;
	if((fr = fopen("residuos.txt","r"))==NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	
	char buf[BUFSIZ];
	int col = 0;
	while(!feof(fr)){
		memset(buf,0,BUFSIZ);
		if(fgets(buf,BUFSIZ,fr)==NULL)break;
		mpz_init_set_str(p, buf, 10);
		int contDiv = 0;
		while((mpz_divisible_p(resultf,p)!=0 && mpz_cmp_ui(p,0)!=0))
		{
			mpz_divexact(resultf,resultf,p);
			contDiv++;
			if(mpz_cmp_ui(resultf,1)==0){
				break;
			}
		}
		data_d[col].col = col;
		data_d[col].n_div = contDiv;
		col++;
		mpz_clear(p);
	}
	fclose(fr);
	if(mpz_cmp_si(resultf,1)==0){
		mpz_clear(resultf);
		return 1;
	}else{
		mpz_clear(resultf);
		return 0;
	}
}


void agregarAVectorDiv(){
	int cont=0;
	while(cont<i_long){		
		insertarNumero(posFila,data_d[cont].col,data_d[cont].n_div%2);
		cont++;
	}
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

