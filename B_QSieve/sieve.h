/** 
* @brief 
* @param n: 
* @param p: numero primo 
* @param r1: root 1 parametro de salida
* @param r2: root 2 parametro de salida
* @return 1 si se encontro solucion o 0 si no se encuentra una solucion
*/
int shanksTonelli(mpz_t n, mpz_t p, mpz_t r1, mpz_t r2);
float *sievingNaive(qs_struct * qs_data, enum TypeSieving typeSieving);
unsigned long *sieving(qs_struct * qs_data, unsigned long *length);
