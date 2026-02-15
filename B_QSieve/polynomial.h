#include <gmp.h>
#include "structsqs.h"

unsigned long fermat(qs_struct * qs_data, unsigned long numLote, unsigned long numPosiciones);
unsigned long standard(qs_struct * qs_data, unsigned long numLote, unsigned long numPosiciones);

// MPQS helpers
int generate_mpqs_poly(qs_struct *qs_data);
int eval_mpqs_Qx(qs_struct *qs_data, mpz_t x, mpz_t result);