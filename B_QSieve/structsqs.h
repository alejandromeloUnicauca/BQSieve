#ifndef STRUCTSQS_H
#define STRUCTSQS_H

/**
 * @file
 * @author Jhon Alejandro Melo<alejandromelo@unicauca.edu.co>
 * @brief Contiene las estructuras para el proyecto
 */
 
 /**
  * @brief estructura que contiene los datos que nos ayudan
  * a recuperar el vector de un numero suave
  */
typedef struct{
	/** maximo comun divisor  */
	mpz_t gcd;
	/** veces que se repite el gcd */
	int periodo;
	/** bloque en el que se obtuvo el gcd y el periodo*/
	int block;
}div_data;

 /**
  * @brief esta estructura contiene un apuntador para crear un
  * array donde se almacenaran los datos que nos ayudan a recuperar el vector
  * exponente de un numero suave
  */
typedef struct{
	/**array de datos de las divisiones*/
	div_data * data;
	/**numero de elementos en el array*/
	int n_values;
}div_data_table;

/**
 * @brief estructura para almacenar los primos
 * de la base y sus logaritmos
 * */
typedef struct{
	mpz_t value;
	mpfr_t log_value;
	unsigned long llog_value;
}prime;

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
} siqs_poly_state;

typedef struct{
    unsigned long large_prime; // primo grande (residuo tras trial division)
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
	/** tabla de parciales para 1LP */
	partials_table partials;
	/** límite para large primes: large_mult * primo_más_grande_de_la_base */
	unsigned long large_prime_bound;
	/** parámetros de criba interpolados para este N */
	sieve_param_t sieve_params;
}qs_struct;

//divisiones sucesivas
typedef struct{
	int col;
	int n_div;
}data_divT;

#endif // STRUCTSQS_H
