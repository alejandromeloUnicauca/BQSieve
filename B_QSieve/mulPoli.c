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
#include <unistd.h>

int main(int argc, char * argv[]){	
	
	mpz_t n;
	mpz_init(n);
	mpz_set_str(n, argv[1], 10);
	printf("N:");
	mpz_out_str(stdout,10,n);
	printf("\n");
	
	FILE * file;
	if ((file = fopen("salida.txt", "r")) == NULL) // open file
	{
		fprintf(stderr,"Falta archivo salida.txt");
		exit(EXIT_FAILURE);
	}
	
	mpz_t qx;
	mpz_t tmp;
	mpz_init(qx);
	mpz_init(tmp);
	
	/*mpz_set_str(qx,"30717132431278347664",10);
	mpz_set_si(tmp,-2818568624543165);
	printf("%d",mpz_congruent_p(qx,tmp,n));
	
	
	exit(EXIT_SUCCESS);*/
	mpz_set_ui(qx,1);
	char buf[BUFSIZ];
	while(!feof(file)){
		memset(buf, 0, BUFSIZ);
		if(fgets(buf,BUFSIZ,file)!=NULL){
			mpz_init_set_str(tmp, buf, 10);
			mpz_mul(qx,qx,tmp);
		}
	}
	fclose(file);
	
	if(mpz_sgn(qx)<0)
		mpz_mul_si(qx,qx,-1);
	/*mpz_out_str(stdout,10,qx);
	printf("\n");*/
	int res = mpz_perfect_square_p(qx);
	printf("%d\n",res);
	
	if(res!=0)
	{
		if ((file = fopen("pos.txt", "r")) == NULL) // open file
		{
			fprintf(stderr,"Falta archivo pos.txt");
			exit(EXIT_FAILURE);
		}
		
		mpz_t mulX;
		mpz_init(mulX);
		mpz_set_ui(mulX,1);
		while(!feof(file))
		{
			memset(buf, 0, BUFSIZ);
			if(fgets(buf,BUFSIZ,file)!=NULL){
				mpz_init_set_str(tmp, buf, 10);
				mpz_mul(mulX,mulX,tmp);
			}
		}
		fclose(file);
		
		/*printf("mul Xi:");
		mpz_out_str(stdout,10,mulX);
		printf("\n\n");*/
		
		mpz_t sqr;
		mpz_init(sqr);
		mpz_t gcd;
		mpz_init(gcd);
		mpz_t res;
		mpz_init(res);
		
		mpz_sqrt(sqr,qx);
		//mpz_pow_ui(mulX,mulX,2);
		
		
		//printf("%d",mpz_perfect_square_p(mulX));
		int con = mpz_congruent_p(qx,mulX,n);
		//printf("con:%d",con);
		/*printf("sqrt:");
		mpz_out_str(stdout,10,sqr);
		printf("\n\n");*/
		
		mpz_sub(res,sqr,mulX);
		
		/*printf("resta:");
		mpz_out_str(stdout,10,res);
		printf("\n\n");*/
		
		mpz_gcd(gcd,res,n);
		
		mpz_t p, q;
		mpz_inits(p,q,NULL);
		
		mpz_set(p,gcd);
		
		/*printf("p:");
		mpz_out_str(stdout,10,gcd);
		printf("\n");*/
		
		mpz_add(res,sqr,mulX);
		
		/*printf("suma:");
		mpz_out_str(stdout,10,res);
		printf("\n\n");*/
		
		mpz_gcd(gcd,res,n);
		
		mpz_set(q,gcd);
		
		/*printf("q:");
		mpz_out_str(stdout,10,gcd);
		printf("\n\n");*/
		
		if(mpz_cmp_ui(p,1)!=0 && mpz_cmp_ui(q,1)!=0){
			gmp_printf("P:%Zd \nQ:%Zd\n",p,q);
			mpz_clears(p,q,NULL);
			exit(EXIT_SUCCESS);
		}
		
		mpz_clear(mulX);
		mpz_clear(sqr);
		mpz_clear(gcd);
		mpz_clear(res);
		
		exit(EXIT_FAILURE);
	}
	
	
	mpz_clear(qx);
	mpz_clear(tmp);
	mpz_clear(n);
	exit(EXIT_FAILURE);
}


