#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <mpfr.h>
#include "sieve.h"
#include "structsqs.h"

int main(int argc, char * argv[]){	
	
	/*mpz_t n, p, r1, r2;
	mpz_inits(n,p,r1,r2,NULL);	
	
	mpz_set_str(n,"15347",10);
	mpz_set_str(p,"31",10);

	shanksTonelli(n,p,r1,r2);
	
	gmp_printf("r1:%Zd\n", r1);
	gmp_printf("r2:%Zd\n", r2);
	fflush(stdout);
	
	mpz_clears(n,p,r1,r2,NULL);*/
	
	qs_struct qs_data;
	
	mpz_inits(qs_data.n,qs_data.interval_length,NULL);
	mpz_set_str(qs_data.n,"15347",10);
	mpz_set_str(qs_data.interval_length,"142",10);
	qs_data.base.length = 5;
	qs_data.base.primes = (prime*)malloc((qs_data.base.length)*sizeof(prime));
	
	mpz_init(qs_data.base.primes[0].value);
	mpz_init(qs_data.base.primes[1].value);
	mpz_init(qs_data.base.primes[2].value);
	mpz_init(qs_data.base.primes[3].value);
	mpz_init(qs_data.base.primes[4].value);
	
	mpz_set_str(qs_data.base.primes[0].value,"2",10);
	mpz_set_str(qs_data.base.primes[1].value,"17",10);
	mpz_set_str(qs_data.base.primes[2].value,"23",10);
	mpz_set_str(qs_data.base.primes[3].value,"29",10);
	mpz_set_str(qs_data.base.primes[4].value,"31",10);
	
	mpfr_t pTemp;
	mpfr_init(pTemp);
	
	mpfr_init(qs_data.base.primes[0].log_value);
	mpfr_init(qs_data.base.primes[1].log_value);
	mpfr_init(qs_data.base.primes[2].log_value);
	mpfr_init(qs_data.base.primes[3].log_value);
	mpfr_init(qs_data.base.primes[4].log_value);
	
	mpfr_set_z(pTemp,qs_data.base.primes[0].value,MPFR_RNDZ);
	mpfr_log(qs_data.base.primes[0].log_value, pTemp, MPFR_RNDZ);
	
	
	mpfr_set_z(pTemp,qs_data.base.primes[1].value,MPFR_RNDZ);
	mpfr_log(qs_data.base.primes[1].log_value, pTemp, MPFR_RNDZ);
	
	mpfr_set_z(pTemp,qs_data.base.primes[2].value,MPFR_RNDZ);
	mpfr_log(qs_data.base.primes[2].log_value, pTemp, MPFR_RNDZ);
	
	mpfr_set_z(pTemp,qs_data.base.primes[3].value,MPFR_RNDZ);
	mpfr_log(qs_data.base.primes[3].log_value, pTemp, MPFR_RNDZ);
	
	mpfr_set_z(pTemp,qs_data.base.primes[4].value,MPFR_RNDZ);
	mpfr_log(qs_data.base.primes[4].log_value, pTemp, MPFR_RNDZ);
	
	
	mpfr_clear(pTemp);
	
	long lengthXi = 0;
	long *Xi = sieving(&qs_data,&lengthXi);
	
	for (long i = 0; i < lengthXi ; i++)
	{
		printf("%ld\n",Xi[i]);
	}
	
	printf("Xis:%ld\n",lengthXi);
	
	mpz_clear(qs_data.base.primes[0].value);
	mpz_clear(qs_data.base.primes[1].value);
	mpz_clear(qs_data.base.primes[2].value);
	mpz_clear(qs_data.base.primes[3].value);
	mpz_clear(qs_data.base.primes[4].value);
	
	mpz_clears(qs_data.n,qs_data.interval_length,NULL);
	exit(EXIT_SUCCESS);
}


