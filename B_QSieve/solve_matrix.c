#include <stdio.h>
#include <stdlib.h>
#include "linalg.h"

int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "Uso: solve_matrix matrix.txt K.sols.txt\n");
        return 1;
    }

    sparse_matrix_t mat;
    if (linalg_read_matrix(argv[1], &mat) != 0)
        return 1;

    int n_words = (mat.n_rows + 63) / 64;
    uint64_t *solutions[LINALG_N_SOLUTIONS];
    for (int k = 0; k < LINALG_N_SOLUTIONS; k++) {
        solutions[k] = calloc(n_words, sizeof(uint64_t));
        if (!solutions[k]) {
            fprintf(stderr, "solve_matrix: sin memoria\n");
            return 1;
        }
    }

    int n_sol = 0;
    linalg_block_lanczos(&mat, solutions, &n_sol);

    int ret = linalg_write_ksols(argv[2], (uint64_t **)solutions, n_sol, mat.n_rows);

    linalg_free_sparse(&mat);
    for (int k = 0; k < LINALG_N_SOLUTIONS; k++)
        free(solutions[k]);

    if (ret != 0) return 1;
    return (n_sol == LINALG_N_SOLUTIONS) ? 0 : 1;
}
