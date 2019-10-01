
/** 
* @brief 
* @param n:
* @param prime:
* @param r1:
* @param r2:
* @return 1 si se encontro solucion o 0 si no se encuentra una solucion
*/
int shanksTonelli(mpz_t n, mpz_t p, mpz_t r1, mpz_t r2);


double *sievingNaivePositive(qs_struct * qs_data);
double *sievingNaiveNegative(qs_struct * qs_data);
long * sieving(qs_struct * qs_data, long *length);
