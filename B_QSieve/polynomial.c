#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "structsqs.h"


/**
 * @brief Calcula los valores de la factorización de Fermat en un rango de posiciones.
 *
 * Esta función calcula los valores de la factorización de Fermat para un rango
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
	unsigned long lastProcessedPosition = 0;
    unsigned long endPosition = numLote * numPosiciones;
    
    if (endPosition > qs_data->intervalo.length_Xi) {
        endPosition = qs_data->intervalo.length_Xi;  // Ajustar si excede el tamaño de los datos
    }else{
		endPosition = numPosiciones;
	}
	
	// Liberar memoria previamente asignada si es necesario
    if (qs_data->intervalo.Qxi != NULL) {
        for (unsigned long i = 0; i < numPosiciones; i++) {
            mpz_clear(qs_data->intervalo.Qxi[i]);
        }
        free(qs_data->intervalo.Qxi);
    }
	
	//Se asigna memoria para el array Qxi
	qs_data->intervalo.Qxi = (mpz_t *)malloc(numPosiciones * sizeof(mpz_t));
	if (qs_data->intervalo.Qxi == NULL) {
        fprintf(stderr, "Error al asignar memoria para Qxi\n");
        exit(EXIT_FAILURE);
    }
    

    for (unsigned long i = 0; i < endPosition; i++) {
		mpz_init(qs_data->intervalo.Qxi[i]);
        //printf("%ld\n,", qs_data->intervalo.Xi[i]);
        mpz_set_si(qs_data->intervalo.Qxi[i], qs_data->intervalo.Xi[posXi]);
        mpz_pow_ui(qs_data->intervalo.Qxi[i], qs_data->intervalo.Qxi[i], 2);
        mpz_sub(qs_data->intervalo.Qxi[i], qs_data->intervalo.Qxi[i], qs_data->n);
        //gmp_printf("%Zd,", qs_data->intervalo.Qxi[i]);
        lastProcessedPosition = i;
        posXi++;
	}
	return lastProcessedPosition;
}

