#include <stdio.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
	//e^(((64/9)^(1/3))*ln(n)^(1/3)*(ln(ln(n)))^(2/3))
	
	mpfr_t num;
	mpfr_init(num);
	mpfr_set_str(num, argv[1], 10, MPFR_RNDU);
	mpfr_t ln1,ln2,e,sqr,pow1,pow2;
	mpfr_init(ln1);
	mpfr_init(ln2);
	mpfr_init(e);
	mpfr_init(sqr);
	mpfr_init(pow1);
	mpfr_init(pow2);
	mpfr_set_str(e, "2.71828182845904523536", 10, MPFR_RNDZ);//define euler
	mpfr_set_str(sqr, "1.9229994270765445097", 10, MPFR_RNDZ);//define  (64/9)^(1/3)
	printf("Raiz:");
	mpfr_out_str(stdout, 10, 0, sqr, MPFR_RNDZ);
	printf("\n:");
	mpfr_set_str(pow1, "0.3333333333333333333", 10, MPFR_RNDZ);//define (1/3)
	mpfr_set_str(pow2, "0.6666666666666666666", 10, MPFR_RNDZ);//define (2/3)
	
	mpfr_log(ln1, num, MPFR_RNDZ); //ln(n)
	printf("ln(n):");
	mpfr_out_str(stdout, 10, 0, ln1, MPFR_RNDZ);
	printf("\n:");
	
	mpfr_log(ln2, ln1, MPFR_RNDZ);//ln(ln(n))
	printf("ln(ln(n)):");
	mpfr_out_str(stdout, 10, 0, ln2, MPFR_RNDZ);
	printf("\n:");
	
	mpfr_pow(ln1, ln1, pow1, MPFR_RNDZ);//ln(n)^(1/3)
	printf("ln(n)^(1/3):");
	mpfr_out_str(stdout, 10, 0, ln1, MPFR_RNDZ);
	printf("\n:");
	
	mpfr_pow(ln2, ln2, pow2, MPFR_RNDZ);//ln(ln(n))^(2/3)
	printf("ln(ln(n))^(2/3):");
	mpfr_out_str(stdout, 10, 0, ln2, MPFR_RNDZ);
	printf("\n:");
	
	mpfr_mul(num, ln1, ln2, MPFR_RNDZ);//ln(n)*(ln(ln(n)))
	printf("ln(ln(n))*ln(n)):");
	mpfr_out_str(stdout, 10, 0, num, MPFR_RNDZ);
	printf("\n:");
	
	mpfr_mul(num, num, sqr, MPFR_RNDZ);//((64/9)^(1/3))*ln(n)*(ln(ln(n)))
	printf("((64/9)^(1/3))*ln(ln(n))*ln(n)):");
	mpfr_out_str(stdout, 10, 0, num, MPFR_RNDZ);
	printf("\n:");
	
	mpfr_pow(num, e, num, MPFR_RNDZ);//e^(((64/9)^(1/3))*ln(n)*(ln(ln(n))))
	
	printf("operaciones:");
	mpfr_out_str(stdout, 10, 0, num, MPFR_RNDZ);
	printf("\n:");
	return 0;
}

