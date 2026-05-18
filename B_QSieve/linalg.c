/*
 * linalg.c — Eliminación Gaussiana GF(2) para la fase de álgebra lineal de MPQS/SIQS.
 *
 * Busca 64 vectores en el nullspace izquierdo de A (n_rows × n_cols) sobre GF(2).
 * Precondición: n_rows = n_cols + 64, rank(A) = n_cols (típico en MPQS).
 *
 * Algoritmo:
 *   1. Convertir A sparse → densa (bitpacked uint64_t).
 *   2. Augmentar con H = I_{n_rows} (tracking de operaciones de fila).
 *   3. Eliminación Gaussiana RREF sobre columnas de A; mismas ops en H.
 *   4. Las n_rows - rank filas restantes de H son vectores del nullspace.
 *
 * Memoria: O(n_rows² / 64) bytes para H.
 *   n_rows=264  (~100 bits): <1 MB   — instantáneo
 *   n_rows=2064 (~183 bits): ~0.5 MB — ~ms
 *   n_rows=27064 (~249 bits): ~110 MB — ~segundos
 * Para N > ~250 bits se necesita Block Lanczos (extensión futura).
 */

#include "linalg.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

#define WORDS(n) (((n) + 63) / 64)

static inline int  bit_get(const uint64_t *row, int j)
    { return (row[j >> 6] >> (j & 63)) & 1; }
static inline void bit_set(uint64_t *row, int j)
    { row[j >> 6] |= (uint64_t)1 << (j & 63); }

static void row_xor(uint64_t *a, const uint64_t *b, int len)
{
    for (int w = 0; w < len; w++) a[w] ^= b[w];
}

/* ======== E/S ============================================================== */

int linalg_read_matrix(const char *filename, sparse_matrix_t *mat)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) { perror(filename); return -1; }
    if (fscanf(fp, "%d %d", &mat->n_rows, &mat->n_cols) != 2) {
        fprintf(stderr, "linalg: cabecera inválida en %s\n", filename);
        fclose(fp); return -1;
    }
    if (mat->n_rows <= mat->n_cols) {
        fprintf(stderr, "linalg: se requiere n_rows > n_cols\n");
        fclose(fp); return -1;
    }
    mat->row_start = malloc((mat->n_rows + 1) * sizeof(int));
    mat->row_len   = malloc(mat->n_rows * sizeof(int));
    int cap = mat->n_rows * 12, total = 0;
    mat->row_data  = malloc(cap * sizeof(int));
    if (!mat->row_start || !mat->row_len || !mat->row_data) {
        fclose(fp); return -1;
    }
    for (int i = 0; i < mat->n_rows; i++) {
        int cnt;
        if (fscanf(fp, "%d", &cnt) != 1) {
            fprintf(stderr, "linalg: error leyendo fila %d\n", i);
            fclose(fp); return -1;
        }
        mat->row_len[i] = cnt;
        mat->row_start[i] = total;
        if (total + cnt > cap) {
            while (total + cnt > cap) cap *= 2;
            mat->row_data = realloc(mat->row_data, cap * sizeof(int));
        }
        for (int k = 0; k < cnt; k++) {
            if (fscanf(fp, "%d", &mat->row_data[total + k]) != 1) {
                fprintf(stderr, "linalg: error leyendo columna\n");
                fclose(fp); return -1;
            }
        }
        total += cnt;
    }
    mat->row_start[mat->n_rows] = total;
    fclose(fp);
    return 0;
}

void linalg_free_sparse(sparse_matrix_t *mat)
{
    free(mat->row_data); free(mat->row_start); free(mat->row_len);
    mat->row_data = mat->row_start = mat->row_len = NULL;
}

/* ======== Solver GF(2) ===================================================== */

int linalg_block_lanczos(const sparse_matrix_t *mat,
                          uint64_t **solutions, int *n_solutions_out)
{
    int nr = mat->n_rows, nc = mat->n_cols;
    int wA = WORDS(nc);
    int wH = WORDS(nr);

    size_t mem_bytes = (size_t)nr * (wA + wH) * sizeof(uint64_t);
    if (mem_bytes >> 20 > 4096) {
        fprintf(stderr, "linalg: matriz demasiado grande (%zu MB). "
                "Se requiere Block Lanczos para N > ~250 bits.\n",
                mem_bytes >> 20);
        *n_solutions_out = 0;
        return -1;
    }

    uint64_t *A = calloc((size_t)nr * wA, sizeof(uint64_t));
    uint64_t *H = calloc((size_t)nr * wH, sizeof(uint64_t));
    if (!A || !H) {
        fprintf(stderr, "linalg: sin memoria (%zu MB)\n", mem_bytes >> 20);
        free(A); free(H);
        *n_solutions_out = 0;
        return -1;
    }

    /* Sparse → densa; H = identidad */
    for (int i = 0; i < nr; i++) {
        bit_set(H + (size_t)i * wH, i);
        int st = mat->row_start[i], ln = mat->row_len[i];
        for (int k = 0; k < ln; k++)
            bit_set(A + (size_t)i * wA, mat->row_data[st + k]);
    }

    /* Eliminación Gaussiana RREF sobre columnas de A */
    int rank = 0;
    for (int col = 0; col < nc && rank < nr; col++) {
        int pivot = -1;
        for (int row = rank; row < nr; row++) {
            if (bit_get(A + (size_t)row * wA, col)) { pivot = row; break; }
        }
        if (pivot < 0) continue;

        if (pivot != rank) {
            uint64_t *rA = A + (size_t)rank  * wA, *pA = A + (size_t)pivot * wA;
            uint64_t *rH = H + (size_t)rank  * wH, *pH = H + (size_t)pivot * wH;
            for (int w = 0; w < wA; w++) { uint64_t t=rA[w]; rA[w]=pA[w]; pA[w]=t; }
            for (int w = 0; w < wH; w++) { uint64_t t=rH[w]; rH[w]=pH[w]; pH[w]=t; }
        }

        uint64_t *pivA = A + (size_t)rank * wA;
        uint64_t *pivH = H + (size_t)rank * wH;
        for (int row = 0; row < nr; row++) {
            if (row == rank) continue;
            if (bit_get(A + (size_t)row * wA, col)) {
                row_xor(A + (size_t)row * wA, pivA, wA);
                row_xor(H + (size_t)row * wH, pivH, wH);
            }
        }
        rank++;
    }

    /* Filas [rank, nr) de H → vectores del nullspace izquierdo */
    int n_sol = 0;
    for (int i = rank; i < nr && n_sol < LINALG_N_SOLUTIONS; i++) {
        uint64_t *hrow = H + (size_t)i * wH;
        uint64_t nz = 0;
        for (int w = 0; w < wH; w++) nz |= hrow[w];
        if (!nz) continue;
        memcpy(solutions[n_sol++], hrow, wH * sizeof(uint64_t));
    }

    *n_solutions_out = n_sol;
    free(A); free(H);

    if (n_sol < LINALG_N_SOLUTIONS)
        fprintf(stderr, "linalg: encontradas %d/%d soluciones\n",
                n_sol, LINALG_N_SOLUTIONS);
    return (n_sol == LINALG_N_SOLUTIONS) ? 0 : -1;
}

/* ======== Escritura de K.sols.txt ========================================== */

int linalg_write_ksols(const char *filename,
                        uint64_t **solutions, int n_solutions, int n_rows)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) { perror(filename); return -1; }
    for (int j = 0; j < n_rows; j++) {
        uint64_t word = 0;
        for (int k = 0; k < n_solutions; k++) {
            if ((solutions[k][j / 64] >> (j % 64)) & 1)
                word |= (uint64_t)1 << (63 - k);
        }
        fprintf(fp, "%016" PRIX64 "\n", word);
    }
    fclose(fp);
    return 0;
}
