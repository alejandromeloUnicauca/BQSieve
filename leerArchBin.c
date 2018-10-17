#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char * argv[]){	
	FILE * fileBin = fopen(argv[1],"rb");
	mpz_t n;//declaracion de variable
	mpz_init(n);//instancia de variable
	char buf[5];
	while(!feof(fileBin)){
		//memset(buf, 1, BUFSIZ);
		fread(buf,1,5,fileBin);
		printf("%s",buf);
		/*mpz_init_set_str(n, buf, 10);
		mpz_out_str(stdout,10,n);*/
	}
	fclose(fileBin);
	return 1;
}


