#ifndef STRUCTSQS_H
#define STRUCTSQS_H

#include <stdint.h>

/**
 * @file
 * @author Jhon Alejandro Melo<alejandromelo@unicauca.edu.co>
 * @brief Contiene las estructuras para el proyecto
 */
 
/**
 * @brief estructura para almacenar los primos
 * de la base y sus logaritmos
 * */
typedef struct{
	mpz_t value;
	mpfr_t log_value;
	unsigned long llog_value;
	uint32_t p;            /* primo como entero nativo */
	uint32_t sqrt_n_mod_p; /* sqrt(N) mod p, precomputado */
}prime;

/**
 * @brief Datos "calientes" de cada primo en un array compacto y contiguo.
 *
 * El bucle de criba y trialDivisionRecip solo necesitan estos campos. Vivir
 * en un array denso de ~20 bytes (en vez de leerlos del struct prime de ~80
 * bytes, que arrastra mpz_t/mpfr_t) mejora mucho la densidad de cache.
 */
typedef struct{
	uint32_t p;        /* primo como entero nativo */
	uint32_t root1;    /* raíz de criba 1 para el polinomio actual */
	uint32_t root2;    /* raíz de criba 2 para el polinomio actual */
	uint32_t recip;    /* ⌊2^32/p⌋ (o +1) para el test de divisibilidad */
	uint8_t  logp;     /* round(log2(p)) */
	uint8_t  rcorrect; /* 1 si recip es exacto, 0 si redondeado arriba */
}sieve_prime;

 /**
  * @brief estructura que hace parte de la tabla de bloques en
  * la que se almacenan los factores y el producto de esos factores
  */
typedef struct{
	/**factores del numero prod_factores*/
	prime * factors;
	/**multiplicacion de todos los factores*/
	mpz_t prod_factors;
	/**cantidad de factores*/
	int length;
}prime_block;

 /**
  * @brief esta estructura contiene un apuntador para crear un
  * array donde se almacenaran los bloques de resiudos cuadraticos
  */
typedef struct{
	/**array de bloques*/
	prime_block * block;
	/**cantidad de bloques*/
	unsigned long length;
}blocks_table;

typedef struct{
	int ** data;
	int n_rows;
	int n_cols;
}matrix;

typedef struct{
	prime * primes;
	sieve_prime * sp;   /* array compacto paralelo a primes[], datos calientes */
	long length;
}primes_base;

typedef struct{
	mpz_t * Xi;
	/**Tamaño del array Xi*/
	unsigned long long length_Xi;
	mpz_t * Qxi;
	/**Tamaño del array Qxi*/
	unsigned long long length_Qxi;
	/**Tamaño del invervalo*/
	mpz_t length;
}interval;

/* Definición de polinomio MPQS colocada antes de qs_struct */
typedef struct{
    mpz_t a;
    mpz_t b;
    mpz_t c;
} mpqs_poly;

/**
 * @brief Máximo número de factores primos que componen 'a' en SIQS.
 */
#define MAX_SIQS_FACTORS 20

/**
 * @brief Estado del generador de polinomios SIQS.
 *
 * Para un a = q_0 * q_1 * ... * q_{s-1}, se generan 2^(s-1) valores
 * de b usando código Gray. Cada b se obtiene sumando/restando 2*B_j
 * al b anterior.
 */
typedef struct {
    int initialized;               /* 0 = necesita nuevo 'a', 1 = iterando sobre b's */
    unsigned int num_factors;      /* s = número de factores de 'a' */
    unsigned long factor_fb_idx[MAX_SIQS_FACTORS]; /* índices en la base de primos */
    mpz_t factors[MAX_SIQS_FACTORS];   /* los primos q_j que componen a */
    mpz_t Bvals[MAX_SIQS_FACTORS];     /* valores auxiliares B_j */
    unsigned long poly_index;      /* índice actual del polinomio derivado (0..2^(s-1)-1) */
    unsigned long num_derived;     /* 2^(s-1) = total de polinomios por este a */
    mpz_t target_a;                /* valor óptimo de a = sqrt(2N)/sieve_size */
    /* Señalización para la actualización Gray code de raíces de criba:
     * sieve_new_a=1 → recomputar raíces y deltas desde cero;
     * sieve_new_a=0 → update incremental con el factor sieve_flip_j y sieve_flip_sign. */
    int sieve_new_a;
    unsigned int sieve_flip_j;
    int sieve_flip_sign;
} siqs_poly_state;

typedef struct{
    unsigned long large_prime;  // primer primo grande (0 si no hay)
    unsigned long large_prime2; // segundo primo grande (0 si 1LP, != 0 si 2LP)
    mpz_t lhs;                // valor a*x+b asociado a esta relación
    mpz_t Qx;                 // valor Q(x) original (con signo)
    mpz_t roota;              // roota del polinomio que generó esta relación
    int *exponents;           // vector de exponentes (tamaño = base.length), sin signo
    int sign;                 // 1 si Q(x) < 0, 0 si Q(x) >= 0
    unsigned int num_a_factors;                  // número de factores de a (SIQS)
    unsigned long a_factor_fb_idx[MAX_SIQS_FACTORS]; // índices en la base de los factores de a
    mpz_t a_value;            // valor de 'a' para esta relación
} partial_entry;

typedef struct{
    partial_entry * entries;
    unsigned long n;
    unsigned long capacity;
} partials_table;

/**
 * @brief Parámetros de criba precompilados, indexados por número de bits de N.
 * Adaptado de la tabla de msieve (Jason Papadopoulos, dominio público).
 */
typedef struct {
	unsigned int bits;       /* tamaño en bits del número a factorizar */
	unsigned int fb_size;    /* tamaño de la base de primos */
	unsigned int large_mult; /* multiplicador para large prime bound (reservado) */
	unsigned int sieve_size; /* mitad del intervalo de criba: se criba [-sieve_size, +sieve_size] */
} sieve_param_t;

/**
 * @brief estructura que contiene los datos necesarios para factorizar
 * un numero n
 * */
typedef struct{
	blocks_table blocks; 
	//TODO:crear estructira para los numeros BSuaves
	long n_BSuaves;
	/***/
	primes_base base;
	/**estrucura para almacenar los vectores exponentes*/
	matrix mat;
	/**Numero que se va factorizar*/
	mpz_t n;
	/**longitud del intervalo positivo*/
	interval intervalo;
	/**polinomio MPQS actual (opcional)*/
	mpqs_poly poly;
	/**raíz roota persistente para generar sucesivos polinomios MPQS (legacy, pre-SIQS)*/
	mpz_t roota;
	/** estado del generador SIQS (múltiples b por cada a) */
	siqs_poly_state siqs_state;
	/** tabla de parciales para 1LP y 2LP */
	partials_table partials;
	/** límite para large primes: large_mult * primo_más_grande_de_la_base */
	unsigned long large_prime_bound;
	/** límite para double large primes: LP_bound^1.8 (cofactor máximo para 2LP) */
	unsigned long long large_prime_bound2;
	/** cuadrado del primo más grande de la base (umbral mínimo para 2LP) */
	unsigned long long max_fb2;
	/** contadores de relaciones 2LP */
	unsigned long n_dlp_stored;
	unsigned long n_dlp_combined;
	/** parámetros de criba interpolados para este N */
	sieve_param_t sieve_params;
	/** multiplicador Knuth-Schroeppel (k tal que factorizamos k*N) */
	unsigned int multiplier;
	/** pool reutilizable para exp_vec — evita calloc por candidato */
	int *exp_vec_pool;
	/** cutoff_config = 1.5 * log2(large_prime_bound), precomputado en main */
	unsigned int sieve_cutoff_config;
	/** pools reutilizables para Xi/Qxi — evita malloc/init/clear por polinomio */
	mpz_t *Xi_pool;
	mpz_t *Qxi_pool;
	unsigned long Xi_Qxi_pool_cap;
	/** hash table O(1) para lookup de parciales 1LP por large prime
	 *  open-addressing con tombstones; keys=0 vacío, ULONG_MAX tombstone */
	unsigned long *lp1_hash_keys;
	unsigned long *lp1_hash_idxs;  /* índice en partials.entries[] */
	unsigned long  lp1_hash_size;  /* potencia de 2 */
	unsigned long  lp1_hash_mask;  /* lp1_hash_size - 1 */
}qs_struct;

#endif // STRUCTSQS_H
