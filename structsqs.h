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
}data_div;

 /**
  * @brief esta estructura contiene un apuntador para crear un
  * array donde se almacenaran los datos que nos ayudan a recuperar el vector
  * exponente de un numero suave
  */
typedef struct{
	/**array de datos de las divisiones*/
	data_div * data;
	/**numero de elementos en el array*/
	int n_values;
}data_div_table;

 /**
  * @brief estructura que hace parte de la tabla de bloques en
  * la que se almacenan los factores y el producto de esos factores
  */
typedef struct{
	/**factores del numero prod_factores*/
	mpz_t * factores;
	/**multiplicacion de todos los factores*/
	mpz_t prod_factores;
	/**cantidad de factores*/
	int n_factores;
}block_remainders;

 /**
  * @brief esta estructura contiene un apuntador para crear un
  * array donde se almacenaran los bloques de resiudos cuadraticos
  */
typedef struct{
	/**array de bloques*/
	block_remainders * blocks;
	/**cantidad de bloques*/
	int n_blocks;
}remainders_block_table;

typedef struct{
	int * data;
	int n_rows;
	int n_cols;
}matrix;

typedef struct{
	int col;
	int n_div;
}data_divT;
