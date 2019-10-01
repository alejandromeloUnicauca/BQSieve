#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <mpfr.h>
#include "structsqs.h"
#include <string.h>

/** 
* @brief 
* @param n: 
* @param p: numero primo 
* @param r1: root 1 parametro de salida
* @param r2: root 2 parametro de salida
* @return 1 si se encontro solucion o 0 si no se encuentra una solucion
*/
int shanksTonelli(mpz_t n, mpz_t p, mpz_t r1, mpz_t r2) {
	
	mpz_t resMod, p1, div;//p1 sera prime-1, resMod y div seran variables temporales para resultados de operaciones
	mpz_inits(resMod,p1,div,NULL);
	
	/*Se puede omitir esta parte por que los primos que llegan 
	 * ya se les calculo el simbolo de legendre
	if((mpz_legendre(n,p)==1)){
		mpz_set_ui(r1,0);
		mpz_set_ui(r2,0);
		return 0;
	}*/
	
	mpz_t q, ss;
	mpz_inits(q,ss,NULL);
	mpz_sub_ui(q,p,1);//q=(prime-1)
	
	
	//mientras que el ultimo bit de q sea 0 (q par)
	while (mpz_divisible_ui_p(q,2) != 0)
	{
		mpz_add_ui(ss,ss,1);
		mpz_divexact_ui(q,q,2);
	}

	if (mpz_cmp_ui(ss,1) == 0)
	{
		mpz_add_ui(p1,p,1);//p1=prime+1;
		mpz_divexact_ui(div,p1,4);
		//div=(prime+1)/4
		mpz_powm(r1,n,div,p);//r1=n^((prime+1)/4) mod prime
		mpz_sub(r2,p,r1);//r2=prime-r1
		
	}else{
		mpz_sub_ui(p1,p,1);//p=prime-1;
		
		mpz_t z;
		mpz_init(z);
		mpz_set_ui(z,2);
		
		mpz_divexact_ui(div,p1,2);
		//(prime-1)/2
		
		mpz_powm(resMod,z,div,p);//resMod=z^(div) mod prime
		
		while (mpz_cmp(resMod,p1)!=0)
		{
			mpz_add_ui(z,z,1);
			mpz_powm(resMod,z,div,p);
		}
		
		mpz_t c,r,t,m;
		mpz_inits(c,r,t,m,NULL);
		
		mpz_powm(c,z,q,p);
		mpz_powm(t,n,q,p);
		
		mpz_add_ui(q,q,1);
		mpz_divexact_ui(div,q,2);//div=q/2
		
		mpz_powm(r,n,div,p);
		
		mpz_set(m,ss);
		
		while (1)
		{
			
			if (mpz_cmp_ui(t,1) == 0)
			{
				mpz_set(r1,r);
				mpz_t pr;
				mpz_init(pr);
				mpz_sub(pr,p,r);
				mpz_set(r2,pr);//r2=prime-r
				mpz_clear(pr);
				break;
			}
			
			mpz_t i,zz,m1;
			mpz_inits(i,zz,m1,NULL);
			
			mpz_set_ui(i,0);
			mpz_set(zz,t);
			
			mpz_sub_ui(m1,m,1);
			
			//zz != 1 && i < (m-1)
			while (mpz_cmp_ui(zz,1) != 0 && mpz_cmp(i,m1) < 0)
			{
				mpz_powm_ui(zz,zz,2,p);//zz=zz*zz mod prime
				mpz_add_ui(i,i,1);
				
			}
			
			mpz_t b,e;
			mpz_inits(b,e,NULL);
			
			mpz_set(b,c);
			mpz_set(e,m);
			mpz_sub(e,e,i);
			mpz_sub_ui(e,e,1);
			
			while (mpz_cmp_ui(e,0) > 0)
			{
				mpz_powm_ui(b,b,2,p);//b=b*b mod prime
				mpz_sub_ui(e,e,1);
			}
			
			mpz_mul(r,r,b);
			mpz_powm_ui(r,r,1,p);
			
			mpz_powm_ui(c,b,2,p);
			
			mpz_mul(t,t,c);
			mpz_powm_ui(t,t,1,p);

			mpz_set(m,i);
		}
	}
	return 1;
}

/** 
* @brief 
* @param qs_data: 
* @return 
*/
double * sievingNaiveNegative(qs_struct * qs_data){
	
	long intervalo = mpz_get_ui(qs_data->interval_length);
	double *S = (double*) malloc(intervalo*sizeof(double));
	memset(S,0,sizeof(double)*intervalo);
	
	mpz_t n,raizn;
	mpz_inits(n,raizn,NULL);
	
	mpz_set(n,qs_data->n);
	mpz_sqrt(raizn,n);
	//gmp_printf("raizn:%Zd",raizn);
	for (int i = 0; i < qs_data->base.length; i++)
	{
		mpz_t x1, x2, p;
		mpz_inits(x1,x2,p,NULL);	
		
		mpz_set(p,qs_data->base.primes[i].value);
		double logp = mpfr_get_d(qs_data->base.primes[i].log_value,MPFR_RNDZ);
		shanksTonelli(n,p,x1,x2);
		//gmp_printf("%Zd,%Zd\n",x1,x2);
		mpz_sub(x1,x1,raizn);
		mpz_mod(x1,x1,p);
		
		mpz_sub(x2,x2,raizn);
		mpz_mod(x2,x2,p);
		//gmp_printf("%Zd,%Zd\n",x1,x2);
		//si x1>x2 se intercambian
		if(mpz_cmp(x1,x2)==1){
			mpz_t temp;
			mpz_init(temp);
			mpz_set(temp,x1);
			mpz_set(x1,x2);
			mpz_set(x2,temp);
			mpz_clear(temp);
		}
		//gmp_printf("%Zd,%Zd\n",x1,x2);
		
		mpz_t d1,d2,Mp;
		mpz_inits(d1,d2,Mp,NULL);
		//x1<x2
		
		mpz_sub(x2,x2,x1);
		
		
		mpz_set(d1,x2);
		mpz_sub(d2,p,d1);

		mpz_sub(Mp,qs_data->interval_length,d1);
		
		//gmp_printf("p:%Zd, x1:%Zd, x2:%Zd, d1:%Zd, d2:%Zd, Mp:%Zd\n",qs_data->base.primes[i].value,x1,x2,d1,d2,Mp);
		
		mpz_t x;
		mpz_init(x);

		mpz_mul_si(Mp,Mp,-1);
		mpz_sub(x1,x1,d2);
		
		//gmp_printf("p:%Zd, x1:%Zd, x2:%Zd, d1:%Zd, d2:%Zd, Mp:%Zd\n",qs_data->base.primes[i].value,x1,x2,d1,d2,Mp);
		for (mpz_set(x,x1); mpz_cmp(x,Mp) == 1 ;)
		{
			S[mpz_get_ui(x)] = S[mpz_get_ui(x)] + logp;
			mpz_sub(x,x,d1);
			S[mpz_get_ui(x)] = S[mpz_get_ui(x)] + logp;
			mpz_sub(x,x,d2);
		}
		
		if(mpz_cmp(x,qs_data->interval_length) == 1){	
			S[mpz_get_ui(x)] = S[mpz_get_ui(x)] + logp;
		}
		
		mpz_clears(x1,x2,x,d1,d2,Mp,p,NULL);
	}
	
	//T=log(sqrt(2N)M)
	
	mpfr_t T;
	mpfr_init_set_ui(T,2,MPFR_RNDZ);
	
	mpfr_mul_z(T,T,qs_data->n,MPFR_RNDZ);
	mpfr_sqrt(T,T,MPFR_RNDZ);
	mpfr_mul_z(T,T,qs_data->interval_length,MPFR_RNDZ);
	mpfr_log(T, T, MPFR_RNDZ);
	mpfr_printf ("T:%.2Rf\n", T);
	
	mpfr_clear(T);
	mpz_clears(n,raizn,NULL);
	return S;
}

/** 
* @brief 
* @param qs_data: 
* @return 
*/
double * sievingNaivePositive(qs_struct * qs_data){
	
	long intervalo = mpz_get_ui(qs_data->interval_length);
	double *S = (double*) malloc(intervalo*sizeof(double));
	memset(S,0,sizeof(double)*intervalo);
	mpz_t n,raizn;
	mpz_inits(n,raizn,NULL);
	
	mpz_set(n,qs_data->n);
	mpz_sqrt(raizn,n);
	//gmp_printf("raizn:%Zd",raizn);
	
	for (int i = 0; i < qs_data->base.length; i++)
	{
		mpz_t x1, x2, p;
		mpz_inits(x1,x2,p,NULL);	
		
		mpz_set(p,qs_data->base.primes[i].value);
		double logp = mpfr_get_d(qs_data->base.primes[i].log_value,MPFR_RNDZ);
		shanksTonelli(n,p,x1,x2);
		//gmp_printf("%Zd,%Zd\n",x1,x2);
		mpz_sub(x1,x1,raizn);
		mpz_mod(x1,x1,p);
		
		mpz_sub(x2,x2,raizn);
		mpz_mod(x2,x2,p);
		//gmp_printf("%Zd,%Zd\n",x1,x2);
		//si x1>x2 se intercambian
		if(mpz_cmp(x1,x2)==1){
			mpz_t temp;
			mpz_init(temp);
			mpz_set(temp,x1);
			mpz_set(x1,x2);
			mpz_set(x2,temp);
			mpz_clear(temp);
		}
		//gmp_printf("%Zd,%Zd\n",x1,x2);
		
		mpz_t d1,d2,Mp;
		mpz_inits(d1,d2,Mp,NULL);
		//x1<x2
		
		mpz_sub(x2,x2,x1);
		
		mpz_set(d1,x2);
		mpz_sub(d2,p,d1);

		mpz_sub(Mp,qs_data->interval_length,d1);
		
		//gmp_printf("p:%Zd, x1:%Zd, x2:%Zd, d1:%Zd, d2:%Zd, Mp:%Zd\n",qs_data->base.primes[i].value,x1,x2,d1,d2,Mp);
		
		mpz_t x;
		mpz_init(x);

		for (mpz_set(x,x1); mpz_cmp(x,Mp) == -1 ;)
		{
			S[mpz_get_ui(x)] = S[mpz_get_ui(x)] + logp;
			mpz_add(x,x,d1);
			S[mpz_get_ui(x)] = S[mpz_get_ui(x)] + logp;
			mpz_add(x,x,d2);
		}
		
		if(mpz_cmp(x,qs_data->interval_length) == -1){	
			S[mpz_get_ui(x)] = S[mpz_get_ui(x)] + logp;
		}
		
		mpz_clears(x1,x2,x,d1,d2,Mp,p,NULL);
	}
	
	//T=log(sqrt(2N)M)
	
	mpfr_t T;
	mpfr_init_set_ui(T,2,MPFR_RNDZ);
	
	mpfr_mul_z(T,T,qs_data->n,MPFR_RNDZ);
	mpfr_sqrt(T,T,MPFR_RNDZ);
	mpfr_mul_z(T,T,qs_data->interval_length,MPFR_RNDZ);
	mpfr_log(T, T, MPFR_RNDZ);
	mpfr_printf ("T:%.2Rf\n", T);
	mpfr_clear(T);
	
	mpz_clears(n,raizn,NULL);
	return S;
}

/** 
* @brief 
* @param qs_data: 
* @param out length:
* @return 
*/
long * sieving(qs_struct * qs_data, long *length){
	
	double *sp = sievingNaivePositive(qs_data);
	double *sn = sievingNaiveNegative(qs_data);
	long intervalLength = mpz_get_ui(qs_data->interval_length);
	long * Xi = (long*)malloc((intervalLength*2)*sizeof(long));
	long contXi=0;
	
	mpz_t raizn;
	mpz_init(raizn);
	mpz_sqrt(raizn,qs_data->n);
	unsigned long raiznl = mpz_get_ui(raizn);
	for (long i = 0; i < intervalLength; i++)
	{
		if(sp[i]>3){
			Xi[contXi] = (i+raiznl);
			contXi++;
			//printf("%ld:%f\n",i,sp[i]);
		}
		
		if(sn[i]>3){
			Xi[contXi] = (-i+raiznl);
			contXi++;
			//printf("-%ld:%f\n",i,sn[i]);
		}
	}
	
	*length = contXi;
	
	free(sp);
	free(sn);
	
	return Xi;
}

	
