#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "structsqs.h"

/* Forward declarations */
int trialDivision(mpz_t Qxi, qs_struct * qs_data, mpz_t Xi);

void insertarNumero(matrix * matriz, int posFila, int posColumna, int valor){

	/*printf("Fila:%d Columna:%d Valor:%d\n",posFila,posColumna,valor);
	printf("Fila:%d Columna:%d Valor:%d\n",posFila,posColumna,matriz->data[posFila][posColumna]);
	fflush(stdout);*/
	
    // Verificar si la matriz y los índices son válidos antes de continuar
    if (matriz == NULL || matriz->data == NULL ||
        posFila < 0 || posFila >= matriz->n_rows ||	
        posColumna < 0 || posColumna >= matriz->n_cols) {
		fprintf(stderr,"Error insertarNumero: Fila=%d (max=%d) Columna=%d (max=%d) Valor=%d\n",
			posFila, matriz->n_rows, posColumna, matriz->n_cols, valor);
		fflush(stderr);
        exit(EXIT_FAILURE);
    }
	
	if(matriz->data != NULL){
		if(matriz->data != NULL && matriz->data[posFila][posColumna] == 0 && valor == 0){
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
}

void agregarAVectorBlock(qs_struct * qs_data, div_data_table * block_table){
	mpz_t gcd;
	mpz_init(gcd);
	for (long i = 0; i < block_table->n_values; i++)
	{
		mpz_set(gcd,block_table->data[i].gcd);
		int block = block_table->data[i].block;
		int periodo = block_table->data[i].periodo;
		
		//si el periodo es par se insertan 0 
		if(periodo & 0){
			for (long j = 0; j < qs_data->blocks.block[block].length; j++)
			{
				insertarNumero(&qs_data->mat,qs_data->n_BSuaves,j+1,0);
			}
		}else{
			mpz_t p;
			mpz_init(p);
			for (long j = 0; j < qs_data->blocks.block[block].length; j++)
			{
				mpz_set(p,qs_data->blocks.block[block].factors[j].value);
				if(mpz_divisible_p(gcd,p)){
					insertarNumero(&qs_data->mat,qs_data->n_BSuaves,(block)*qs_data->blocks.block[0].length+j+1,periodo%2);
				}else{
					insertarNumero(&qs_data->mat,qs_data->n_BSuaves,(block)*qs_data->blocks.block[0].length+j+1,0);
				}
			}
			mpz_clear(p);
		}
	}
	mpz_clear(gcd);
}

int blockDivision(mpz_t Qxi, qs_struct * qs_data){
	div_data_table block_table;
	block_table.data = (div_data*)malloc((qs_data->base.length+1)*sizeof(div_data));
	mpz_t QxiTemp;
	mpz_init(QxiTemp);
	mpz_set(QxiTemp,Qxi);
	if(mpz_sgn(QxiTemp)==-1)
		mpz_mul_si(QxiTemp,QxiTemp,-1);
		
	mpz_t gcd, gcdAnt;
	mpz_inits(gcd,gcdAnt,NULL);
	unsigned long contGcd = 0, cont = 0;
	for (long i = 0; i < qs_data->blocks.length; i++)
	{
		contGcd = 0;
		mpz_gcd(gcd,QxiTemp,qs_data->blocks.block[i].prod_factors);
		mpz_set(gcdAnt,gcd);
		if (mpz_cmp_ui(gcd,1) != 0) {
			while(mpz_cmp_ui(gcd,1)!=0){
				mpz_divexact(QxiTemp,QxiTemp,gcd);
				contGcd++;
				mpz_gcd(gcd,QxiTemp,qs_data->blocks.block[i].prod_factors);
				//si el maximo comun divisor cambia se guardan los datos en la tabla
				if(mpz_cmp(gcd,gcdAnt)!=0){
					//Se almacena en una estructura el gcd,
					//las veces que se repite, y el bloque al que pertenece
					block_table.data[cont].block = i;
					mpz_init(block_table.data[cont].gcd);
					mpz_set(block_table.data[cont].gcd,gcdAnt);
					block_table.data[cont].periodo = contGcd;
					mpz_set(gcdAnt,gcd);
					contGcd = 0;
					cont++;
				}
			}
		}
	}

	block_table.n_values = cont;
	if(mpz_cmp_ui(QxiTemp,1)==0){
		agregarAVectorBlock(qs_data, &block_table);
		mpz_clears(QxiTemp,gcd,gcdAnt,NULL);
		for (unsigned long k = 0; k < cont; k++){
			mpz_clear(block_table.data[k].gcd);
		}
		free(block_table.data);
		return 1;
	}else{
		mpz_clears(QxiTemp,gcd,gcdAnt,NULL);
		for (unsigned long k = 0; k < cont; k++){
			mpz_clear(block_table.data[k].gcd);
		}
		free(block_table.data);
		return 0;
	}
}

/**
 * @brief Factoriza el array Qxi con bloques. Cuando blockDivision falla,
 * intenta trialDivision como fallback para capturar large primes.
 * @param qs_data estructura que contiene el array Qxi
 * @param endPos cantidad de candidatos
 * @param posXi índice inicial
 * @return retorna 1 si aun faltan numeros B_suaves por verificar y 0 en caso de haberlos encontrado todos
 */
int factoringBlocks(qs_struct * qs_data,  unsigned long endPos, unsigned long posXi){

	FILE * fp;
	if((fp = fopen("polinomio.txt","a")) == NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	
	for (unsigned long i = 0; i < endPos; i++)
	{
		if(blockDivision(qs_data->intervalo.Qxi[i],qs_data)==1){
			/* Full relation via bloques */
			if(mpz_sgn(qs_data->intervalo.Qxi[i]) < 0){
				insertarNumero(&qs_data->mat, qs_data->n_BSuaves, 0, 1);
			}
			qs_data->n_BSuaves++;
			mpz_t lhs;
			mpz_init(lhs);
			mpz_mul(lhs, qs_data->poly.a, qs_data->intervalo.Xi[posXi]);
			mpz_add(lhs, lhs, qs_data->poly.b);
			mpz_out_str(fp, 10, lhs);
			fprintf(fp, ";");
			mpz_out_str(fp, 10, qs_data->intervalo.Qxi[i]);
			fprintf(fp, ";");
			mpz_out_str(fp, 10, qs_data->roota);
			fprintf(fp, "\n");
			fflush(fp);
			mpz_clear(lhs);
			if(qs_data->n_BSuaves==qs_data->mat.n_rows){
				fclose(fp);
				return 0;
			}
		} else {
			/* blockDivision falló: intentar trialDivision para capturar 1LP */
			int result = trialDivision(qs_data->intervalo.Qxi[i], qs_data, qs_data->intervalo.Xi[posXi]);
			if(result == 1){
				/* Full relation via trial (raro aquí, pero posible) */
				qs_data->n_BSuaves++;
				mpz_t lhs;
				mpz_init(lhs);
				mpz_mul(lhs, qs_data->poly.a, qs_data->intervalo.Xi[posXi]);
				mpz_add(lhs, lhs, qs_data->poly.b);
				mpz_out_str(fp, 10, lhs);
				fprintf(fp, ";");
				mpz_out_str(fp, 10, qs_data->intervalo.Qxi[i]);
				fprintf(fp, ";");
				mpz_out_str(fp, 10, qs_data->roota);
				fprintf(fp, "\n");
				fflush(fp);
				mpz_clear(lhs);
				if(qs_data->n_BSuaves==qs_data->mat.n_rows){
					fclose(fp);
					return 0;
				}
			} else if(result == 2){
				/* Combined partial — ya escrito por trialDivision */
				qs_data->n_BSuaves++;
				if(qs_data->n_BSuaves==qs_data->mat.n_rows){
					fclose(fp);
					return 0;
				}
			}
			/* result == 0: parcial guardada o descartada */
		}
		posXi++;
	}
	fclose(fp); 
	return 1;
}

void agregarAVectorDiv(qs_struct * qs_data, data_divT * data_d){

	for (long i = 0; i < qs_data->base.length ; i++)
	{
		insertarNumero(&qs_data->mat,qs_data->n_BSuaves,data_d[i].col+1,data_d[i].n_div%2);
	}
}



/**
 * @brief Valida si un numero del polinomio se divide en la base de residuos
 * usando divisiones triviales. Si el residuo es 1, es B-suave (retorna 1).
 * Si el residuo es un primo < large_prime_bound, se guarda como parcial
 * y se intenta emparejar (retorna 2 si se combinó, 0 si solo se guardó).
 * @param Qxi: Numero que se valida si es divisible en la base
 * @param qs_data: estructura con base de primos y tabla de parciales
 * @param Xi: valor x del candidato (para calcular lhs = a*x+b)
 * @return 1 = full relation, 2 = combined partial, 0 = no relation
 */
int trialDivision(mpz_t Qxi, qs_struct * qs_data, mpz_t Xi){
	int *exp_vec = (int*)calloc(qs_data->base.length, sizeof(int));
	unsigned long contDiv = 0;
	mpz_t QxiTemp;
	mpz_init(QxiTemp);
	mpz_set(QxiTemp,Qxi);
	int sign = 0;
	if(mpz_sgn(QxiTemp)==-1){
		mpz_mul_si(QxiTemp,QxiTemp,-1);
		sign = 1;
	}

	for (unsigned long i = 0; i < qs_data->base.length; i++){
		mpz_t p;
		mpz_init(p);
		mpz_set(p,qs_data->base.primes[i].value);
		contDiv = 0;

		while(mpz_divisible_p(QxiTemp,p)!=0){
			mpz_divexact(QxiTemp,QxiTemp,p);
			contDiv++;
			if(mpz_cmp_ui(QxiTemp,1)==0) break;
		}
		exp_vec[i] = contDiv;
		mpz_clear(p);
	}
	
	if(mpz_cmp_si(QxiTemp,1)==0){
		/* Full relation: insertar vector en la matriz */
		if(sign)
			insertarNumero(&qs_data->mat, qs_data->n_BSuaves, 0, 1);
		for (long i = 0; i < qs_data->base.length; i++)
			insertarNumero(&qs_data->mat, qs_data->n_BSuaves, i+1, exp_vec[i] % 2);
		mpz_clear(QxiTemp);
		free(exp_vec);
		return 1;
	}
	
	/* ¿Es un posible large prime? */
	unsigned long residuo = 0;
	if(mpz_fits_ulong_p(QxiTemp))
		residuo = mpz_get_ui(QxiTemp);
	
	if(residuo > 1 && residuo < qs_data->large_prime_bound &&
	   mpz_probab_prime_p(QxiTemp, 15) > 0) {
		/* Buscar si ya tenemos una parcial con el mismo large prime */
		long match_idx = -1;
		for(unsigned long k = 0; k < qs_data->partials.n; k++){
			if(qs_data->partials.entries[k].large_prime == residuo){
				match_idx = (long)k;
				break;
			}
		}
		
		if(match_idx >= 0){
			/* ¡Match! Combinar las dos parciales para crear una full relation.
			 * Si parcial_1 tiene Q1(x1) = (-1)^s1 * prod(pi^ei) * LP
			 * y parcial_2 tiene Q2(x2) = (-1)^s2 * prod(pi^fi) * LP
			 * entonces Q1*Q2 = (-1)^(s1+s2) * prod(pi^(ei+fi)) * LP²
			 * El vector de exponentes mod 2 es XOR (suma mod 2) de ambos vectores.
			 * LP² es par, así que LP desaparece de la paridad.
			 */
			partial_entry *match = &qs_data->partials.entries[match_idx];
			
			/* Insertar vector combinado (XOR) en la matriz */
			int combined_sign = (sign + match->sign) % 2;
			if(combined_sign)
				insertarNumero(&qs_data->mat, qs_data->n_BSuaves, 0, 1);
			for(long i = 0; i < qs_data->base.length; i++){
				int combined_exp = (exp_vec[i] + match->exponents[i]) % 2;
				insertarNumero(&qs_data->mat, qs_data->n_BSuaves, i+1, combined_exp);
			}
			
			/* Escribir relación combinada en polinomio.txt:
			 * lhs = lhs1 * lhs2 (mod N), Qx = Q1 * Q2
			 * Necesitamos guardar ambos lhs y ambos Q para la raíz cuadrada.
			 * Formato: lhs1*lhs2;Q1*Q2;roota1,roota2
			 * mulPoli multiplicará todos los Q's y todos los lhs's.
			 */
			FILE *fp = fopen("polinomio.txt", "a");
			if(fp){
				/* lhs combinado = lhs1 * lhs2 */
				mpz_t combined_lhs, combined_Q;
				mpz_inits(combined_lhs, combined_Q, NULL);
				
				/* lhs del candidato actual */
				mpz_t my_lhs;
				mpz_init(my_lhs);
				mpz_mul(my_lhs, qs_data->poly.a, Xi);
				mpz_add(my_lhs, my_lhs, qs_data->poly.b);
				
				/* lhs combinado */
				mpz_mul(combined_lhs, my_lhs, match->lhs);
				
				/* Q combinado = Q1 * Q2 (NO dividimos por LP²;
				 * LP² tiene exponente par y no afecta la paridad
				 * en la matriz. mulPoli lo incluirá en el producto
				 * y será un cuadrado perfecto junto con lhs1*lhs2) */
				mpz_mul(combined_Q, Qxi, match->Qx);
				
				/* Escribir: lhs_combinado;Q_combinado;roota_actual,roota_match */
				mpz_out_str(fp, 10, combined_lhs);
				fprintf(fp, ";");
				mpz_out_str(fp, 10, combined_Q);
				fprintf(fp, ";");
				mpz_out_str(fp, 10, qs_data->roota);
				fprintf(fp, ",");
				mpz_out_str(fp, 10, match->roota);
				fprintf(fp, "\n");
				fflush(fp);
				fclose(fp);
				
				mpz_clears(combined_lhs, combined_Q, my_lhs, NULL);
			}
			
			/* Eliminar la parcial usada (swap con última) */
			free(match->exponents);
			mpz_clears(match->lhs, match->Qx, match->roota, NULL);
			unsigned long last = qs_data->partials.n - 1;
			if((unsigned long)match_idx != last)
				qs_data->partials.entries[match_idx] = qs_data->partials.entries[last];
			qs_data->partials.n--;
			
			mpz_clear(QxiTemp);
			free(exp_vec);
			return 2; /* combined partial = full relation */
		} else {
			/* No hay match: guardar esta parcial */
			if(qs_data->partials.n >= qs_data->partials.capacity){
				unsigned long newcap = qs_data->partials.capacity == 0 ? 1024 : qs_data->partials.capacity * 2;
				qs_data->partials.entries = realloc(qs_data->partials.entries, newcap * sizeof(partial_entry));
				qs_data->partials.capacity = newcap;
			}
			unsigned long idx = qs_data->partials.n++;
			partial_entry *e = &qs_data->partials.entries[idx];
			e->large_prime = residuo;
			e->exponents = exp_vec; /* transferir ownership */
			e->sign = sign;
			mpz_init(e->lhs);
			mpz_mul(e->lhs, qs_data->poly.a, Xi);
			mpz_add(e->lhs, e->lhs, qs_data->poly.b);
			mpz_init_set(e->Qx, Qxi);
			mpz_init_set(e->roota, qs_data->roota);
			
			mpz_clear(QxiTemp);
			/* NO free exp_vec, se transfirió a la parcial */
			return 0;
		}
	}
	
	mpz_clear(QxiTemp);
	free(exp_vec);
	return 0;
}

/**
 * @brief Factoriza el array Qxi con divisiones triviales, verifica si cada posicion es un numero bsuave
 * y lo agrega al archivo polinomio.txt
 * @param qs_data estructura que contiene el array Qxi
 * @param endPos cantidad de candidatos
 * @param posXi índice inicial
 * @return retorna 1 si aun faltan numeros B_suaves por verificar y 0 en caso de haberlos encontrado todos
 */
int factoringTrial(qs_struct * qs_data, unsigned long endPos, unsigned long posXi){
	FILE * fp;
	if((fp = fopen("polinomio.txt","a")) == NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	for (unsigned long i = 0; i < endPos; i++)
	{
		int result = trialDivision(qs_data->intervalo.Qxi[i], qs_data, qs_data->intervalo.Xi[posXi]);
		if(result == 1){
			/* Full relation encontrada — vector ya insertado por trialDivision */
			qs_data->n_BSuaves++;
			/* Escribir (a*x+b);Q(x);roota */
			mpz_t lhs;
			mpz_init(lhs);
			mpz_mul(lhs, qs_data->poly.a, qs_data->intervalo.Xi[posXi]);
			mpz_add(lhs, lhs, qs_data->poly.b);
			mpz_out_str(fp, 10, lhs);
			fprintf(fp, ";");
			mpz_out_str(fp, 10, qs_data->intervalo.Qxi[i]);
			fprintf(fp, ";");
			mpz_out_str(fp, 10, qs_data->roota);
			fprintf(fp, "\n");
			fflush(fp);
			mpz_clear(lhs);
			if(qs_data->n_BSuaves==qs_data->mat.n_rows){
				fclose(fp);
				return 0;
			}
		} else if(result == 2){
			/* Combined partial — vector y polinomio.txt ya escritos por trialDivision */
			qs_data->n_BSuaves++;
			if(qs_data->n_BSuaves==qs_data->mat.n_rows){
				fclose(fp);
				return 0;
			}
		}
		/* result == 0: no relation o parcial guardada */
		posXi++;
	}
	fclose(fp); 
	return 1;
}


