void insertarNumero(matrix * matriz, int posFila, int posColumna, int valor);
void add_a_factors_to_matrix(qs_struct *qs_data);
int blockDivisionV2(mpz_t Qxi, qs_struct *qs_data, mpz_t Xi);
int trialDivisionRecip(mpz_t Qxi, qs_struct *qs_data, mpz_t Xi, unsigned long sieve_offset);
int factoringTrial(qs_struct * qs_data, unsigned long endPos, unsigned long posXi, unsigned long xmax);
int factoringBlocks(qs_struct * qs_data, unsigned long endPos, unsigned long posXi, unsigned long xmax);

/* Almacén en memoria de relaciones: reemplaza polinomio.txt. */
typedef struct {
    mpz_t lhs_mpz;       /* a*x+b como mpz_t                          */
    mpz_t qfile_mpz;     /* Q(x) = lhs^2 - kN como mpz_t             */
    mpz_t roota_mpz[2];  /* rootas: [0] siempre válido, [1] solo 2LP */
    int   n_roota;       /* número de rootas inicializadas (1 o 2)   */
} rel_entry_t;

typedef struct {
    rel_entry_t *entries;
    int          count;
    int          cap;
} relation_store_t;

void              relstore_init(void);
void              relstore_free(void);
relation_store_t *relstore_get(void);

