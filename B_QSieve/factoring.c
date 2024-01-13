#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "structsqs.h"


void insertarNumero(matrix * matriz, int posFila, int posColumna, int valor){

	/*printf("Fila:%d Columna:%d Valor:%d\n",posFila,posColumna,valor);
	printf("Fila:%d Columna:%d Valor:%d\n",posFila,posColumna,matriz->data[posFila][posColumna]);
	fflush(stdout);*/
	
    // Verificar si la matriz y los índices son válidos antes de continuar
    if (matriz == NULL || matriz->data == NULL ||
        posFila < 0 || posFila >= matriz->n_rows ||	
        posColumna < 0 || posColumna >= matriz->n_cols) {
		//printf("Fila:%d Columna:%d Valor:%d\n",posFila,posColumna,valor);
		//printf("Fila:%d Columna:%d Valor:%d\n",posFila,posColumna,matriz->data[posFila][posColumna]);
		fflush(stdout);
        printf("Error: Parámetros no válidos en insertarNumero\n");
        exit(EXIT_FAILURE);
    }
	
	if(matriz->data != NULL){
		if(matriz->data != NULL && matriz->data[posFila][posColumna] == 0 && valor == 0){
			matriz->data[posFila][posColumna] = 0;
			return;
		}
		
		if(matriz->data[posFila][posColumna] == 1 && valor == 1){
			matriz->data[posFila][posColumna] = 0;
			return;
		}
		
		if(matriz->data[posFila][posColumna] == 0 && valor == 1){
			matriz->data[posFila][posColumna] = 1;
			return;
		}
		
		if(matriz->data[posFila][posColumna] == 1 && valor == 0){
			matriz->data[posFila][posColumna] = 1;
			return;
		}
	}
}

void agregarAVectorBlock(qs_struct * qs_data, div_data_table * block_table){
	mpz_t gcd;
	mpz_init(gcd);
	for (long i = 0; i < block_table->n_values; i++)
	{
		mpz_set(gcd,block_table->data[i].gcd);
		int block = block_table->data[i].block;
		int periodo = block_table->data[i].periodo;
		
		//si el periodo es par se insertan 0 
		if(periodo & 0){
			for (long j = 0; j < qs_data->blocks.block[block].length; j++)
			{
				insertarNumero(&qs_data->mat,qs_data->n_BSuaves,j,0);
			}
		}else{
			mpz_t p;
			mpz_init(p);
			for (long j = 0; j < qs_data->blocks.block[block].length; j++)
			{
				mpz_set(p,qs_data->blocks.block[block].factors[j].value);
				if(mpz_divisible_p(gcd,p)){
					insertarNumero(&qs_data->mat,qs_data->n_BSuaves,(block)*qs_data->blocks.block[0].length+j,periodo%2);
				}else{
					insertarNumero(&qs_data->mat,qs_data->n_BSuaves,(block)*qs_data->blocks.block[0].length+j,0);
				}
			}
			mpz_clear(p);
		}
	}
	mpz_clear(gcd);
}

int blockDivision(mpz_t Qxi, qs_struct * qs_data){
	//TODO:optimizar memoria
	div_data_table block_table;
	block_table.data = (div_data*)malloc((qs_data->base.length+1)*sizeof(div_data));
	mpz_t QxiTemp;
	mpz_init(QxiTemp);
	mpz_set(QxiTemp,Qxi);
	if(mpz_sgn(QxiTemp)==-1)
		mpz_mul_si(QxiTemp,QxiTemp,-1);
		
	mpz_t gcd, gcdAnt;
	mpz_inits(gcd,gcdAnt,NULL);
	long contGcd = 0, cont = 0;
	for (long i = 0; i < qs_data->blocks.length; i++)
	{
		contGcd = 0;
		mpz_gcd(gcd,QxiTemp,qs_data->blocks.block[i].prod_factors);
		mpz_set(gcdAnt,gcd);
		
		while(mpz_cmp_ui(gcd,1)!=0){
			mpz_divexact(QxiTemp,QxiTemp,gcd);
			contGcd++;
			mpz_gcd(gcd,QxiTemp,qs_data->blocks.block[i].prod_factors);
			//si el maximo comun divisor cambia se guardan los datos en la tabla
			if(mpz_cmp(gcd,gcdAnt)!=0){
				//Se almacena en una estructura el gcd,
				//las veces que se repite, y el bloque al que pertenece
				block_table.data[cont].block = i;
				mpz_init(block_table.data[cont].gcd);
				mpz_set(block_table.data[cont].gcd,gcdAnt);
				block_table.data[cont].periodo = contGcd;
				mpz_set(gcdAnt,gcd);
				contGcd = 0;
				cont++;
			}
		}
	}

	block_table.n_values = cont;
	if(mpz_cmp_ui(QxiTemp,1)==0){
		agregarAVectorBlock(qs_data, &block_table);
		mpz_clears(QxiTemp,gcd,gcdAnt,NULL);
		//Todo Liberar memoria de data->gcd
		free(block_table.data);
		return 1;
	}else{
		free(block_table.data);
		mpz_clears(QxiTemp,gcd,gcdAnt,NULL);
		return 0;
	}
}

/**
 * @brief Factoriza el array Qxi con divisiones triviales, verifica si cada posicion es un numero bsuave
 * y lo agrega al archivo polinomio.txt
 * @param qs_data estructura que contiene el array Qxi
 * @param startPos Posición donde inicia la verificación de los números bsuaves
 * @param endPos Posición donde finaliza la verificacion de los números bsuaves
 * @return retorna 1 si aun faltan numeros B_suaves por verificar y 0 en caso de haberlos encontrado todos
 */
int factoringBlocks(qs_struct * qs_data,  unsigned long endPos, unsigned long posXi){

	FILE * fp;//file residuos
	if((fp = fopen("polinomio.txt","a")) == NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	
	for (unsigned long i = 0; i < endPos; i++)
	{
		if(blockDivision(qs_data->intervalo.Qxi[i],qs_data)==1){	
			qs_data->n_BSuaves++;
			if(qs_data->n_BSuaves==qs_data->base.length+1)return 0;
			fprintf(fp,"%ld;",qs_data->intervalo.Xi[posXi]); 
			mpz_out_str(fp,10,qs_data->intervalo.Qxi[i]);
			fprintf(fp,"\n");
		}
		posXi++;
	}
	fclose(fp); 
	return 1;
}

void agregarAVectorDiv(qs_struct * qs_data, data_divT * data_d){

	for (long i = 0; i < qs_data->base.length ; i++)
	{
		insertarNumero(&qs_data->mat,qs_data->n_BSuaves,data_d[i].col,data_d[i].n_div%2);
	}
}



/**
 * @brief Valida si un numero del polinomio se divide en la base de residuos
 * usando el metodo de los bloques, si el numero es liso se almacenan sus factores en un 
 * archivo
 * @param Qxi: Numero que se valida si es divisible en la base
 * @param qs_data: base de primos
 * @return retorna 1 si es divisible en la base o 0 si no lo es
 */
int trialDivision(mpz_t Qxi, qs_struct * qs_data){
	data_divT * data_d;
	data_d = (data_divT*)malloc((qs_data->base.length)*sizeof(data_divT));
	int contDiv = 0;
	mpz_t QxiTemp;
	mpz_init(QxiTemp);
	mpz_set(QxiTemp,Qxi);
	if(mpz_sgn(QxiTemp)==-1)
		mpz_mul_si(QxiTemp,QxiTemp,-1);

	for (unsigned long i = 0; i < qs_data->base.length;)
	{
		mpz_t p;
		mpz_init(p);
		mpz_set(p,qs_data->base.primes[i].value);
		contDiv = 0;

		while(mpz_divisible_p(QxiTemp,p)!=0){
			mpz_divexact(QxiTemp,QxiTemp,p);
			contDiv++;
			if(mpz_cmp_ui(QxiTemp,1)==0){
				break;
			}
		}
		data_d[i].col = i;
		data_d[i].n_div = contDiv;
		i++;
		mpz_clear(p);
	}
	
	if(mpz_cmp_si(QxiTemp,1)==0){
		mpz_clear(QxiTemp);
		agregarAVectorDiv(qs_data,data_d);
		free(data_d);
		return 1;
	}else{
		mpz_clear(QxiTemp);
		free(data_d);
		return 0;
	}
}

/**
 * @brief Factoriza el array Qxi con divisiones triviales, verifica si cada posicion es un numero bsuave
 * y lo agrega al archivo polinomio.txt
 * @param qs_data estructura que contiene el array Qxi
 * @param startPos Posición donde inicia la verificación de los números bsuaves
 * @param endPos Posición donde finaliza la verificacion de los números bsuaves
 * @return retorna 1 si aun faltan numeros B_suaves por verificar y 0 en caso de haberlos encontrado todos
 */
int factoringTrial(qs_struct * qs_data, unsigned long endPos, unsigned long posXi){
	FILE * fp;//file residuos
	if((fp = fopen("polinomio.txt","a")) == NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	for (unsigned long i = 0; i < endPos; i++)
	{
		if(trialDivision(qs_data->intervalo.Qxi[i],qs_data)==1){
			qs_data->n_BSuaves++;
			if(qs_data->n_BSuaves==qs_data->base.length+1)return 0;
			fprintf(fp,"%ld;",qs_data->intervalo.Xi[posXi]); 
			mpz_out_str(fp,10,qs_data->intervalo.Qxi[i]);
			fprintf(fp,"\n");
		}
		posXi++;
	}
	fclose(fp); 
	return 1;
}


