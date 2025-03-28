#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "structsqs.h"


/**
 * @brief Calcula los valores de la factorización de Fermat en un rango de posiciones.
 *
 * Esta función calcula los valores de la factorización de Fermat (x_i^2-N) para un rango
 * de posiciones especificado, comenzando desde la posición de inicio y procesando
 * la cantidad de posiciones indicada. Los resultados se almacenan en el arreglo Qxi
 * dentro de la estructura qs_data.
 * @param qs_data     Estructura que contiene los datos necesarios para el cálculo.
 * @param startPos    Posición de inicio desde donde comenzar a calcular.
 * @param endPos	  Ultima Posicion de posiciones a procesar.
 *
 * @return La última posición calculada dentro del rango especificado.
 */
unsigned long fermat(qs_struct *qs_data, unsigned long numLote, unsigned long numPosiciones){
	unsigned long posXi = (numLote-1)*numPosiciones;
    unsigned long endPosition = numLote * numPosiciones;
    // Liberar memoria previamente asignada si es necesario
    if (qs_data->intervalo.Qxi != NULL) {
        for (unsigned long i = 0; i < numPosiciones; i++) {
            mpz_clear(qs_data->intervalo.Qxi[i]);
        }

        free(qs_data->intervalo.Qxi);
        qs_data->intervalo.Qxi = NULL;
    }

    if (endPosition > qs_data->intervalo.length_Xi) {
        numPosiciones = qs_data->intervalo.length_Xi - posXi ;  // Ajustar si excede el tamaño de los datos
        endPosition = numPosiciones;
    }
    
	//Se asigna memoria para el array Qxi
	qs_data->intervalo.Qxi = (mpz_t *)malloc(numPosiciones * sizeof(mpz_t));
    qs_data->intervalo.length_Qxi = numPosiciones;
	if (qs_data->intervalo.Qxi == NULL) {
        fprintf(stderr, "Error al asignar memoria para Qxi\n");
        exit(EXIT_FAILURE);
    }
    
    unsigned long i;

    for (i = 0; i < numPosiciones; i++) {
		mpz_init(qs_data->intervalo.Qxi[i]);
        mpz_set(qs_data->intervalo.Qxi[i], qs_data->intervalo.Xi[posXi]);
        mpz_pow_ui(qs_data->intervalo.Qxi[i], qs_data->intervalo.Qxi[i], 2);
        mpz_sub(qs_data->intervalo.Qxi[i],qs_data->intervalo.Qxi[i],  qs_data->n);
        // gmp_printf("%Zd\n", qs_data->intervalo.Xi[posXi]);
        posXi++;
	}

	return posXi;
}

/**
 * @brief Calcula los valores de la factorización de Fermat en un rango de posiciones.
 *
 * Esta función calcula los valores de la factorización de Fermat ((x+sqrt(⌈N⌉)^2)-N) para un rango
 * de posiciones especificado, comenzando desde la posición de inicio y procesando
 * la cantidad de posiciones indicada. Los resultados se almacenan en el arreglo Qxi
 * dentro de la estructura qs_data.
 * @param qs_data     Estructura que contiene los datos necesarios para el cálculo.
 * @param startPos    Posición de inicio desde donde comenzar a calcular.
 * @param endPos	  Ultima Posicion de posiciones a procesar.
 *
 * @return La última posición calculada dentro del rango especificado.
 */
unsigned long standard(qs_struct *qs_data, unsigned long numLote, unsigned long numPosiciones){
    unsigned long posXi = (numLote-1)*numPosiciones;
    unsigned long endPosition = numLote * numPosiciones;
    
    // Liberar memoria previamente asignada si es necesario
    if (qs_data->intervalo.Qxi != NULL) {
        for (unsigned long i = 0; i < numPosiciones; i++) {
            mpz_clear(qs_data->intervalo.Qxi[i]);
        }

        free(qs_data->intervalo.Qxi);
        qs_data->intervalo.Qxi = NULL;
    }

    if (endPosition > qs_data->intervalo.length_Xi) {
        numPosiciones = endPosition = qs_data->intervalo.length_Xi;  // Ajustar si excede el tamaño de los datos
    }
	
	//Se asigna memoria para el array Qxi
	qs_data->intervalo.Qxi = (mpz_t *)malloc(numPosiciones * sizeof(mpz_t));
    qs_data->intervalo.length_Qxi = numPosiciones;
	if (qs_data->intervalo.Qxi == NULL) {
        fprintf(stderr, "Error al asignar memoria para Qxi\n");
        exit(EXIT_FAILURE);
    }
    
    unsigned long i;

    for (i = 0; i < numPosiciones; i++) {
        mpz_t c_sqrtN;
        mpfr_t sqrtN;
        mpz_init(c_sqrtN);
		mpz_init(qs_data->intervalo.Qxi[i]);
        mpfr_init2(sqrtN,mpz_sizeinbase(qs_data->n,2));
        mpfr_set_z(sqrtN,qs_data->n,MPFR_RNDN);
        mpz_set(qs_data->intervalo.Qxi[i], qs_data->intervalo.Xi[posXi]);
        mpfr_sqrt(sqrtN,sqrtN,MPFR_RNDZ);
        mpfr_get_z(c_sqrtN,sqrtN,MPFR_RNDU);
        mpz_add(qs_data->intervalo.Qxi[i], qs_data->intervalo.Qxi[i], c_sqrtN);
        mpz_pow_ui(qs_data->intervalo.Qxi[i],qs_data->intervalo.Qxi[i],2);
        mpz_sub(qs_data->intervalo.Qxi[i],qs_data->intervalo.Qxi[i],qs_data->n);
        posXi++;
        mpz_clear(c_sqrtN);
        mpfr_clear(sqrtN);
	}

	return endPosition;
}

