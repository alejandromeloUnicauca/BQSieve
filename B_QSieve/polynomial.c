#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "structsqs.h"


void fermat(qs_struct * qs_data){
	//TODO:mejorar manejo de memoria
	mpz_t * Qxi;	
	Qxi = (mpz_t*)malloc((qs_data->intervalo.length_Xi)*sizeof(mpz_t));
	
	for (long i = 0; i < qs_data->intervalo.length_Xi; i++)
	{
		mpz_init(Qxi[i]);
		//printf("%ld,",qs_data->intervalo.Xi[i]);
		mpz_set_si(Qxi[i],qs_data->intervalo.Xi[i]);
		mpz_pow_ui(Qxi[i],Qxi[i],2);
		mpz_sub(Qxi[i],Qxi[i],qs_data->n);
		//gmp_printf("%Zd,",Qxi[i]);
	}
	
	qs_data->intervalo.Qxi = Qxi;
	qs_data->intervalo.length_Qxi = qs_data->intervalo.length_Xi;
}

