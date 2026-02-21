/** 
* @brief 
* @param n: 
* @param p: numero primo 
* @param r1: root 1 parametro de salida
* @param r2: root 2 parametro de salida
* @return 1 si se encontro solucion o 0 si no se encuentra una solucion
*/
int shanksTonelli(mpz_t n, mpz_t p, mpz_t r1, mpz_t r2);

/**
 * @brief Precomputa sqrt(N) mod p y campos nativos para cada primo de la base.
 * Se llama UNA VEZ tras generar la base de primos.
 */
void sieve_precompute_roots(qs_struct *qs_data);

/**
 * @brief Calcula raíces de criba para el polinomio actual.
 * Se llama cada vez que se genera un nuevo polinomio.
 */
void sieve_compute_roots(qs_struct *qs_data, unsigned long sieve_interval);

/**
 * @brief Criba logarítmica uint8 en bloques de 32KB (cache-friendly).
 * Calcula las raíces de criba para Q(x)=a*x²+2*b*x+c y acumula log(p)
 * en bloques de uint8. Devuelve los índices que superan el umbral.
 *
 * @param qs_data   Estructura con la base de primos y polinomio actual
 * @param xmax      Mitad del intervalo de criba (se criba [-xmax..xmax])
 * @param out_indices  Array de salida con las posiciones x que pasan la criba
 * @param out_count    Número de candidatos encontrados
 */
void sieve_mpqs(qs_struct *qs_data, unsigned long xmax,
                long **out_indices, unsigned long *out_count);
