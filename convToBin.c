#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char * argv[]){	
	FILE * fileASC = fopen(argv[1],"r");
	FILE * fileBin = fopen(argv[2],"wb");
	mpz_t num;//declaracion de variable
	mpz_init(num);//instancia de variable
	char buf[BUFSIZ];
	while(!feof(fileASC)){
		memset(buf, 0, BUFSIZ);
		fgets(buf,BUFSIZ,fileASC);
		printf("%s",buf);
		mpz_init_set_str(num, buf, 10);
		//int nchar = strlen(buf);
		//mpz_out_str(fileBin,2,&num);
		fwrite(&buf,sizeof(char),strlen(buf),fileBin);
		//fwrite(&num,sizeof(mpz_t),1,fileBin);
	}
	fclose(fileASC);
	fclose(fileBin);
	
	return 1;
}

