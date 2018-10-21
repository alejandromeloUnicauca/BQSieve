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
	int gcd;
	/** veces que se repite el gcd */
	int periodo;
	/** bloque en el que se obtuvo el gcd y el periodo*/
	int block;
}blocks_div_table;
