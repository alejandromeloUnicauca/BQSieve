void insertarNumero(matrix * matriz, int posFila, int posColumna, int valor);
void add_a_factors_to_matrix(qs_struct *qs_data);
int blockDivisionV2(mpz_t Qxi, qs_struct *qs_data, mpz_t Xi);
int trialDivisionRecip(mpz_t Qxi, qs_struct *qs_data, mpz_t Xi, unsigned long sieve_offset);
int factoringTrial(qs_struct * qs_data, unsigned long endPos, unsigned long posXi, unsigned long xmax);
int factoringBlocks(qs_struct * qs_data, unsigned long endPos, unsigned long posXi, unsigned long xmax);

/* Almacén en memoria de relaciones: reemplaza polinomio.txt. */
typedef struct {
    char *lhs;    /* cadena decimal de a*x+b                       */
    char *qfile;  /* cadena decimal de Q(x) = lhs^2 - kN           */
    char *rootas; /* rootas separadas por coma, o NULL si ninguna   */
} rel_entry_t;

typedef struct {
    rel_entry_t *entries;
    int          count;
    int          cap;
} relation_store_t;

void              relstore_init(void);
void              relstore_free(void);
relation_store_t *relstore_get(void);

