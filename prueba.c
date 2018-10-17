#include <stdio.h>
#include <gmp.h>

int main(){
  mpz_t n;
  mpz_t n2;
  mpz_init(n);
  mpz_init(n2);
  mpz_set_ui(n,3);
  mpz_set_ui(n2,10);
  if(mpz_divisible_p(n2,n)){
    printf("true");
  }else{
    printf("false");
  }
}
