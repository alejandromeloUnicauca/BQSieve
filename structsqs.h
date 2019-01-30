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
	long gcd;
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
  * @brief estructura que hace parte de la tabla de bloques en
  * la que se almacenan los factores y el producto de esos factores
  */
typedef struct{
	/**factores del numero prod_factores*/
	mpz_t * factors;
	/**multiplicacion de todos los factores*/
	mpz_t prod_factors;
	/**cantidad de factores*/
	int n_factors;
}remainders_block;

 /**
  * @brief esta estructura contiene un apuntador para crear un
  * array donde se almacenaran los bloques de resiudos cuadraticos
  */
typedef struct{
	/**array de bloques*/
	remainders_block * blocks;
	/**cantidad de bloques*/
	int n_blocks;
}blocks_table;

typedef struct{
	int ** data;
	int n_rows;
	int n_cols;
}matrix;

typedef struct{
	blocks_table blocks; 
	mpz_t n_remainders;
	div_data_table blocks_data;
	matrix mat;
	mpz_t n;
	mpz_t interval;
}qs_struct;

//divisiones sucesivas
typedef struct{
	int col;
	int n_div;
}data_divT;
