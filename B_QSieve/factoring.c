#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "structsqs.h"


void insertarNumero(matrix * matriz, int posFila, int posColumna, int valor){

	if(matriz->data[posFila][posColumna] == 0 && valor == 0){
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
	for (long i = 0; i < qs_data->base.length;)
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

void factoringTrial(qs_struct * qs_data){
	FILE * fp;//file residuos
	if((fp = fopen("polinomio.txt","w")) == NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
	}

	for (unsigned long i = 0; i < qs_data->intervalo.length_Qxi; i++)
	{
		if(trialDivision(qs_data->intervalo.Qxi[i],qs_data)==1){
			qs_data->n_BSuaves++;
			if(qs_data->n_BSuaves==qs_data->base.length+1)break;
			fprintf(fp,"%ld;",qs_data->intervalo.Xi[i]); 
			mpz_out_str(fp,10,qs_data->intervalo.Qxi[i]);
			fprintf(fp,"\n");
		}
	}
	fclose(fp); 
}


