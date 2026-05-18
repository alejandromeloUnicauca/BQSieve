#ifndef LINALG_H
#define LINALG_H

#include <stdint.h>
#include <stdio.h>

/*
 * Número de soluciones que produce el solver.
 * Debe coincidir con el excedente de filas sobre columnas en matrix.txt:
 *   n_rows - n_cols = 64  (construido así en B_QSieve.c: n_rows = base.length+1+64)
 */
#define LINALG_N_SOLUTIONS 64

/*
 * sparse_matrix_t — Matriz binaria GF(2) en formato sparse de entrada.
 *
 * row_data: pool plano contiguo de índices de columnas (todos los 1s de todas
 *           las filas concatenados). row_start[i] es el índice en row_data
 *           donde comienzan los 1s de la fila i; row_start[n_rows] apunta al
 *           final del pool (facilita iterar sin row_len).
 * row_len:  número de 1s en la fila i (= row_start[i+1] - row_start[i]).
 */
typedef struct {
    int  *row_data;    /* pool: indices de columnas con bit=1 */
    int  *row_start;   /* row_start[0..n_rows], tamaño n_rows+1 */
    int  *row_len;     /* row_len[0..n_rows-1]                  */
    int   n_rows;      /* = base.length + 1 + 64                */
    int   n_cols;      /* = base.length + 1                     */
} sparse_matrix_t;

/*
 * linalg_read_matrix — Lee matrix.txt en formato sparse ASCII.
 *
 * Formato esperado:
 *   Línea 1: "n_rows n_cols"
 *   Líneas 2..n_rows+1: "count col1 col2 ... colCount"
 *
 * Retorna 0 en éxito, -1 en error (imprime mensaje a stderr).
 */
int linalg_read_matrix(const char *filename, sparse_matrix_t *mat);

/*
 * linalg_free_sparse — Libera la memoria de una sparse_matrix_t.
 */
void linalg_free_sparse(sparse_matrix_t *mat);

/*
 * linalg_block_lanczos — Block Lanczos sobre GF(2) con paralelismo OpenMP.
 *
 * Encuentra el nullspace izquierdo de mat: vectores y ∈ F_2^{n_rows} tal que
 * y^T * A = 0 (mod 2), equivalente a A^T * y = 0.
 *
 * Algoritmo: Montgomery Block Lanczos (1995) con bloque de 64 bits sobre la
 * matriz simétrica B = A × A^T. Itera V_{k+1} = B×V_k + V_k×D_k + V_{k-1}×E_k
 * comprobando en cada paso si alguna columna del bloque cae en null(A^T).
 * Las multiplicaciones matriz-vector dispersas se paralelizan con OpenMP.
 *
 * Complejidad: O(n_cols/64 × nnz) tiempo; O(n_rows × 8) bytes de memoria.
 *   n=10,000 : ~0.01 s, ~1 MB
 *   n=50,000 : ~0.3 s, ~6 MB
 *   n=140,000: ~2 s, ~18 MB
 *
 * solutions[k] debe ser un array pre-allocado de ceil(n_rows/64) uint64_t.
 * Escribe en solutions[0..n_sol-1] los vectores nulos encontrados.
 * Escribe en *n_solutions_out el número de soluciones (esperado: 64).
 *
 * Retorna 0 en éxito, -1 si no se encontraron suficientes null vectors.
 */
int linalg_block_lanczos(const sparse_matrix_t *mat,
                          uint64_t **solutions,
                          int *n_solutions_out);

/*
 * linalg_write_ksols — Escribe K.sols.txt en el formato que consume B_Sieve.py.
 *
 * Produce n_rows líneas, cada una un hex de 16 dígitos (64 bits).
 * El bit (63-k) de la línea j está activo si la relación j pertenece a la
 * solución k. Este mapeo hace que process_polynomial() en B_Sieve.py acceda
 * a la solución k en la iteración i=k+1 (bin_num_padded[k] = bit 63-k).
 *
 * Retorna 0 en éxito, -1 en error de escritura.
 */
int linalg_write_ksols(const char *filename,
                       uint64_t **solutions,
                       int n_solutions,
                       int n_rows);

#endif /* LINALG_H */
