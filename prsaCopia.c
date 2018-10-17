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

void longitudBase(mpfr_t n, mpz_t result);
void intervaloPolinomio(mpfr_t n, mpz_t result);
void calcularResiduos(mpz_t n, mpz_t longitud, mpz_t residuos);
void polinomioFermat(mpz_t n, mpz_t intervalo);
void polinomioFermatB(mpz_t n, mpz_t intervalo);
int validarResultado(mpz_t result);
int validarResultadoB(mpz_t result);
void bloques(int n);
void reducirPolinomio();


int main(int argc, char * argv[]){	
	
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
	mpz_init(intervaloP);
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
	
	//intervalo del polinomio
	mpfr_set_str(num, argv[1], 10, MPFR_RNDU);
	intervaloPolinomio(num,intervaloP);
	printf("intervalo del polinomio: ");
	mpz_out_str(stdout,10,intervaloP);
	printf("\n");
	
	t_inicio = clock();
	//generar residuos
	printf("Generando residuos cuadraticos...\n");
	mpz_t residuos;
	calcularResiduos(n,longitudB,residuos);
	mpz_out_str(stdout,10,residuos);
	printf(" residuos generados en residuos.txt\n");
	t_final = clock();
	segundos = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
	printf("tiempo de calculo de residuos:%fs\n",segundos);
	
	//crear los bloques
	bloques(3);
	
	//Calcular polinomio de Fermat
	t_inicio = clock();
	printf("Calculando polinomio de Fermat\n");
	polinomioFermat(n,intervaloP);
	t_final = clock();
	segundos = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
	printf("tiempo de calculo del polinomio:%fs\n",segundos);
	
	//Calcular polinomio de Fermat por metodo de bloques
	t_inicio = clock();
	printf("Calculando polinomio de Fermat con reduccion por bloques\n");
	polinomioFermatB(n,intervaloP);
	t_final = clock();
	segundos = (double) (t_final-t_inicio)/CLOCKS_PER_SEC;
	printf("tiempo de calculo del polinomio:%fs\n",segundos);


	
	//Liberar espacio de variables
	mpz_clear(n);
	mpz_clear(intervaloP);
	mpz_clear(longitudB);
	mpfr_clear(num);

	exit(EXIT_SUCCESS);
}


/**
 * @brief calcula los residuos cuadraticos del numero n que se obtienen
 * del archivo primes.txt y los almacena en un archivo residuos.txt 
 * @param n: numero que se le evaluaran los residuos cuadraticos
 * @param longitud: numero de residuos que se requieren para el numero n
 * @param residuos: contador que se utiliza para devolver el numero de 
 * residuos encontrados
 */
void calcularResiduos(mpz_t n, mpz_t longitud, mpz_t residuos){
	//result = n^exp % mod
	mpz_init(residuos);//contador de residuos encontrados
	mpz_set_ui(residuos,0);//residuos = 0
	
	mpz_t result, mod, exp;
	
	mpz_init(result);
	mpz_init(mod);
	mpz_init(exp);
	
	FILE * file;//file primes
	FILE * fr;//file residuos
	if((fr = fopen("residuos.txt","w")) == NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	fclose(fr);
	
	mpz_t contn;//contador de primos que se ultilizan de primes.txt
	mpz_init(contn);
	mpz_set_ui(contn,0);
	
	//si el archivo primes.txt no existe termina
	if ((file = fopen("primes.txt", "r")) == NULL) // open file
	{
		fprintf(stderr,"Falta archivo primes.txt");
		exit(EXIT_FAILURE);
	}
	
	char buf[BUFSIZ];
	while(!feof(file)){
		mpz_add_ui(contn,contn,1);//contn = contn + 1
		memset(buf, 0, BUFSIZ);
		if(fgets(buf,BUFSIZ,file)!=NULL){
			mpz_init_set_str(mod, buf, 10);//se asigna un numero del archivo a mod
			mpz_init_set_str(exp, buf, 10);//se asigna un numero del archivo a exp
			mpz_sub_ui(exp,exp,1);//exp = exp-1
			mpz_divexact_ui(exp,exp,2);//exp = exp/2
			mpz_powm(result,n,exp,mod);//potenciacion modular result = (n^exp)%mod
			//si el resultado es 1 se agrega al archivo
			if((mpz_cmp_si(result,1)==0)){
				mpz_add_ui(residuos,residuos,1);//residuos = residuos + 1
				if((fr = fopen("residuos.txt","a"))==NULL){
					perror("fopen");
					exit(EXIT_FAILURE);
				}
				fprintf(fr,"%s",buf);
				fclose(fr);
			}
			//si los residuos es igual a la longitud termina
			if((mpz_cmp(residuos, longitud)==0))
				break;
		}
	}
	
	mpz_clear(result);
	mpz_clear(exp);
	mpz_clear(mod);
	printf("Se usaron ");  
	mpz_out_str(stdout,10,contn);
	printf(" primos\n");
	mpz_clear(contn);
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
 * @brief Calcula los valores del polinomio y los agrega al archivo polinomio.txt
 * @param n: Numero al que se le calculan los valores del polinomio
 * @param intervalo: Nos indica el rango en el que deben evaluarse los residuos del numero n
 */
void polinomioFermat(mpz_t n, mpz_t intervalo){
	//Q(Xi)=(sqrt(n)+i)^2-n
	mpz_t raiz;//sqrt(n)
	mpz_t limite;//limite del ciclo
	mpz_t result;//resultados del polinomio
	
	mpz_init(raiz);
	mpz_init(limite);
	mpz_init(result);
	
	mpz_set(limite,intervalo);//limite = intervalo
	mpz_sqrt(raiz,n);//raiz = sqrt(n)
	mpz_mul_si(intervalo,intervalo,-1);//intervalo = -1*intervalo
	
	FILE * fp;
	if((fp = fopen("polinomio.txt","w"))==NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	fclose(fp);

	char buf[BUFSIZ];
	
	while(mpz_cmp(intervalo,limite)!=1){
		//calcular resultado
		mpz_add(result,raiz,intervalo);//result=raiz+intervalo
		mpz_pow_ui(result,result,2);//result = result^2
		mpz_sub(result,result,n);//result = result-n
		
		//Si al validar el resultado da 1 se agrega al archivo polinomio
		if(validarResultado(result)){
			//agregar al archivo
			memset(buf,0,BUFSIZ);
			if((fp = fopen("polinomio.txt","a"))==NULL){
				perror("fopen");
				exit(EXIT_FAILURE);
			}
			mpz_out_str(fp,10,result);
			fprintf(fp,"\n");
			fclose(fp);
		}
		//se incrementa el valor del intervalo en 1
		mpz_add_ui(intervalo,intervalo,1);
	}
	
	mpz_clear(raiz);
	mpz_clear(limite);
	mpz_clear(result);
}

/**
 * @brief Calcula los valores del polinomio y los agrega al archivo polinomio.txt con el metodo de los bloques
 * @param n: Numero al que se le calculan los valores del polinomio
 * @param intervalo: Nos indica el rango en el que deben evaluarse los residuos del numero n
 */
void polinomioFermatB(mpz_t n, mpz_t intervalo){
	//Q(Xi)=(sqrt(n)+i)^2-n
	mpz_t raiz;//sqrt(n)
	mpz_t limite;//limite del ciclo
	mpz_t result;//resultados del polinomio
	
	mpz_init(raiz);
	mpz_init(limite);
	mpz_init(result);
	
	mpz_set(limite,intervalo);//limite = intervalo
	mpz_sqrt(raiz,n);//raiz = sqrt(n)
	mpz_mul_si(intervalo,intervalo,-1);//intervalo = -1*intervalo
	
	FILE * fp;
	if((fp = fopen("polinomioR.txt","w"))==NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	fclose(fp);

	char buf[BUFSIZ];
	
	while(mpz_cmp(intervalo,limite)!=1){
		//calcular resultado
		mpz_add(result,raiz,intervalo);//result=raiz+intervalo
		mpz_pow_ui(result,result,2);//result = result^2
		mpz_sub(result,result,n);//result = result-n
		
		//Si al validar el resultado da 1 se agrega al archivo polinomio
		if(validarResultadoB(result)){
			//agregar al archivo
			memset(buf,0,BUFSIZ);
			if((fp = fopen("polinomioR.txt","a"))==NULL){
				perror("fopen");
				exit(EXIT_FAILURE);
			}
			mpz_out_str(fp,10,result);
			fprintf(fp,"\n");
			fclose(fp);
		}
		//se incrementa el valor del intervalo en 1
		mpz_add_ui(intervalo,intervalo,1);
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
	mpz_t p;//primo
	mpz_t resultf;
	
	mpz_init(resultf);//variable de resultado final
	mpz_init(p);
	
	mpz_set(resultf,result);//resultf = result
	//si el resultado es negativo se vuelve positivo
	if(mpz_sgn(resultf)==-1)
		mpz_mul_si(resultf,resultf,-1);
	FILE * fr;
	if((fr = fopen("residuos.txt","r"))==NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	char buf[BUFSIZ];
	while(!feof(fr)){
		memset(buf,0,BUFSIZ);
		fgets(buf,BUFSIZ,fr);
		mpz_init_set_str(p, buf, 10);
		while((mpz_divisible_p(resultf,p)!=0 && mpz_cmp_ui(p,0)!=0))
		{
			mpz_divexact(resultf,resultf,p);
			if(mpz_cmp_ui(resultf,1)==0){
				break;
			}
		}
	}
	fclose(fr);
	if(mpz_cmp_si(resultf,1)==0){
		return 1;
	}else{
		return 0;
	}
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
int validarResultadoB(mpz_t result){

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
	
	//TODO:Cambiar archivo temporal por demora
	FILE * fgcTemp;
	fgcTemp = fopen("gcdTemp.txt","w");
	fclose(fgcTemp);
	FILE * fb;
	fb = fopen("bloques.txt","r");
	while(!feof(fb)){
		memset(buf,0,BUFSIZ);
		
		if(fgets(buf,BUFSIZ,fb)==NULL){
			break;
		}
		
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
				fgcTemp = fopen("gcdTemp.txt","a");
				fprintf(fgcTemp,"numero:");
				mpz_out_str(fgcTemp,10,result);
				fprintf(fgcTemp,", mcd:");
				mpz_out_str(fgcTemp,10,mcdA);		
				fprintf(fgcTemp,";%d\n",contGcda);
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

/*void reducirPolinomio(){
	FILE * fp;
	fp = fopen("polinomio.txt","r");

	FILE * fpr;
	fpr = fopen("polinomioR.txt","w");
	fclose(fpr);

	char buf[BUFSIZ];

	mpz_t valorP;//valor del polinomio
	mpz_t valorPC;//copia del valor del polinomio
	mpz_init(valorP);
	mpz_init(valorPC);

	mpz_t valorB;//valor del bloque
	mpz_init(valorB);

	mpz_t mcd;//valor del bloque
	mpz_init(mcd);

	while(!feof(fp)){
		memset(buf,0,BUFSIZ);
		if(fgets(buf,BUFSIZ,fp)!=NULL){
			mpz_init_set_str(valorP, buf, 10);
			mpz_init_set_str(valorPC, buf, 10);
			FILE * fb;
			fb = fopen("bloques.txt","r");
			while(!feof(fb)){
				memset(buf,0,BUFSIZ);
				if(fgets(buf,BUFSIZ,fb)!=NULL){
					mpz_init_set_str(valorB, buf, 10);
					do{//TODO: ARREGLAR CONDICION
						printf("MCD:");
						mpz_gcd(mcd,valorP,valorB);
						mpz_out_str(stdout,10,mcd);
						printf("\n");
						if(mpz_cmp_ui(mcd,1)==0){
							printf("break\n");
							break;
						}
						mpz_divexact(valorP,valorP,mcd);
						printf("Resultado Division:");
						mpz_out_str(stdout,10,valorP);
						printf("\n");
					}while(mpz_cmp_ui(mcd,1)!=0);
				}
			}
			printf("Resultado Final:");
			mpz_out_str(stdout,10,valorP);
			printf("\n");
			fclose(fb);
			if(mpz_cmp_ui(valorP,1)==0){
				fpr = fopen("polinomioR.txt","a");
				mpz_out_str(fpr,10,valorPC);
				fprintf(fpr,"\n");
				fclose(fpr);
			}
		}
	}

	fclose(fp);
}*/
