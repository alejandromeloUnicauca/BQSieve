void insertarNumero(matrix * matriz, int posFila, int posColumna, int valor);
void agregarAVectorBlock(qs_struct * qs_data, div_data_table * block_table);
void agregarAVectorDiv(qs_struct * qs_data, data_divT * data_d);
int blockDivision(mpz_t Qxi, qs_struct * qs_data);
int trialDivision(mpz_t Qxi, qs_struct * qs_data);
void factoringTrial(qs_struct * qs_data);
void factoringBlocks(qs_struct * qs_data);

