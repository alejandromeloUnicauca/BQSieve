#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>

int shanksTonelli(mpz_t n, mpz_t prime, mpz_t r1, mpz_t r2);
void BSNaive();


int main(int argc, char * argv[]){	
	
	mpz_t n, prime, r1, r2;
	mpz_inits(n,prime,r1,r2,NULL);	
	
	mpz_set_str(n,"15347",10);
	mpz_set_str(prime,"31",10);

	shanksTonelli(n,prime,r1,r2);
	
	gmp_printf("r1:%Zd\n", r1);
	gmp_printf("r2:%Zd\n", r2);
	fflush(stdout);
	
	mpz_clears(n,prime,r1,r2,NULL);
	exit(EXIT_SUCCESS);
}

/** 
* @brief 
* @param n:
* @param prime:
* @param r1:
* @param r2:
* @return 1 si se encontro solucion o 0 si no se encuentra una solucion
*/
int shanksTonelli(mpz_t n, mpz_t prime, mpz_t r1, mpz_t r2) {
	
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
	mpz_sub_ui(q,prime,1);//q=(prime-1)
	
	
	//mientras que el ultimo bit de q sea 0 (q par)
	while (mpz_divisible_ui_p(q,2) != 0)
	{
		mpz_add_ui(ss,ss,1);
		mpz_divexact_ui(q,q,2);
	}

	if (mpz_cmp_ui(ss,1) == 0)
	{
		mpz_add_ui(p1,prime,1);//p1=prime+1;
		mpz_divexact_ui(div,p1,4);
		//div=(prime+1)/4
		mpz_powm(r1,n,div,prime);//r1=n^((prime+1)/4) mod prime
		mpz_sub(r2,prime,r1);//r2=prime-r1
		
	}else{
		mpz_sub_ui(p1,prime,1);//p=prime-1;
		
		mpz_t z;
		mpz_init(z);
		mpz_set_ui(z,2);
		
		mpz_divexact_ui(div,p1,2);
		//(prime-1)/2
		
		mpz_powm(resMod,z,div,prime);//resMod=z^(div) mod prime
		
		while (mpz_cmp(resMod,p1)!=0)
		{
			mpz_add_ui(z,z,1);
			mpz_powm(resMod,z,div,prime);
		}
		
		mpz_t c,r,t,m;
		mpz_inits(c,r,t,m,NULL);
		
		mpz_powm(c,z,q,prime);
		mpz_powm(t,n,q,prime);
		
		mpz_add_ui(q,q,1);
		mpz_divexact_ui(div,q,2);//div=q/2
		
		mpz_powm(r,n,div,prime);
		
		mpz_set(m,ss);
		
		while (1)
		{
			
			if (mpz_cmp_ui(t,1) == 0)
			{
				mpz_set(r1,r);
				mpz_t pr;
				mpz_init(pr);
				mpz_sub(pr,prime,r);
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
				mpz_powm_ui(zz,zz,2,prime);//zz=zz*zz mod prime
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
				mpz_powm_ui(b,b,2,prime);//b=b*b mod prime
				mpz_sub_ui(e,e,1);
			}
			
			mpz_mul(r,r,b);
			mpz_powm_ui(r,r,1,prime);
			
			mpz_powm_ui(c,b,2,prime);
			
			mpz_mul(t,t,c);
			mpz_powm_ui(t,t,1,prime);

			mpz_set(m,i);
		}
	}
	return 1;
}

