void insertarNumero(matrix * matriz, int posFila, int posColumna, int valor);
void agregarAVectorBlock(qs_struct * qs_data, div_data_table * block_table);
void agregarAVectorDiv(qs_struct * qs_data, data_divT * data_d);
void add_a_factors_to_matrix(qs_struct *qs_data);
int blockDivision(mpz_t Qxi, qs_struct * qs_data);
int trialDivision(mpz_t Qxi, qs_struct * qs_data, mpz_t Xi);
int trialDivisionRecip(mpz_t Qxi, qs_struct *qs_data, mpz_t Xi, unsigned long sieve_offset);
int factoringTrial(qs_struct * qs_data, unsigned long endPos, unsigned long posXi, unsigned long xmax);
int factoringBlocks(qs_struct * qs_data, unsigned long endPos, unsigned long posXi, unsigned long xmax);

