/**
 * @file
 * @author Jhon Alejandro Melo <alejandromelo@unicauca.edu.co>
 * 			Juan Manuel Campo <>
 * @copyright GNU Public License.
 *
 * @brief Fase de raíz cuadrada del Quadratic Sieve.
 *
 * ==============================================================================
 *  TODO: REFACTOR ESTRUCTURAL PENDIENTE (opción 3 — estilo msieve)
 * ==============================================================================
 *  Actualmente esta fase construye el producto ∏Q_i con miles de bits y le saca
 *  raíz cuadrada. Aún con producto en árbol balanceado (O(N log N)), el
 *  mpz_sqrt sobre un número de >>100 000 bits sigue creciendo y limita la
 *  escalabilidad para claves >256 bits.
 *
 *  Solución correcta: NUNCA construir ∏Q. En su lugar:
 *    1. Durante la criba (factoring.c), persistir el VECTOR DE EXPONENTES
 *       ENTERO de cada relación (no solo el mod 2 que va a la matriz).
 *       Para combined partials, guardar la suma de los vectores enteros.
 *    2. Cambiar formato de polinomio.txt: añadir columna con exp_vec enteros
 *       (o archivo paralelo polinomio.exp).
 *    3. En mulPoli: sumar los exp_vec enteros del subconjunto seleccionado
 *       → vector E con todas las entradas pares (porque ∏Q es cuadrado).
 *       Reconstruir Y = ∏ p_i^(E_i/2) mod N usando mpz_powm_ui — nunca se
 *       construye un número > |N|.
 *    4. mulX = ∏ lhs_i mod N (ya implementado abajo).
 *    5. gcd(Y − mulX, N) y gcd(Y + mulX, N) como ahora.
 *
 *  Esfuerzo: ~150 LOC repartido entre B_QSieve.c, factoring.c, mulPoli.c.
 *  Speedup esperado: 50-100× sobre la versión actual; ESENCIAL para 256+ bits.
 * ==============================================================================
 */

#include <stdio.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void load_mpz_file(const char *path, mpz_t **arr, size_t *count, size_t *cap) {
    FILE *fp = fopen(path, "r");
    if (fp == NULL) return;
    char buf[BUFSIZ];
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        if (buf[0] == '\n' || buf[0] == '\0') continue;
        if (*count == *cap) {
            *cap = (*cap == 0) ? 64 : (*cap * 2);
            *arr = (mpz_t *)realloc(*arr, (*cap) * sizeof(mpz_t));
        }
        mpz_init_set_str((*arr)[*count], buf, 10);
        (*count)++;
    }
    fclose(fp);
}

/* Lee roota_list.txt: cada línea trae uno o varios rootas separados por coma
 * (combined partials). Añade cada roota dos veces porque qx se multiplica por
 * roota². */
static void load_rootas_squared(const char *path, mpz_t **arr, size_t *count, size_t *cap) {
    FILE *fp = fopen(path, "r");
    if (fp == NULL) return;
    char buf[BUFSIZ];
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        if (buf[0] == '\n' || buf[0] == '\0') continue;
        char *token = strtok(buf, ",\n");
        while (token != NULL) {
            while (*token == ' ') token++;
            if (*token != '\0') {
                for (int rep = 0; rep < 2; rep++) {
                    if (*count == *cap) {
                        *cap = (*cap == 0) ? 64 : (*cap * 2);
                        *arr = (mpz_t *)realloc(*arr, (*cap) * sizeof(mpz_t));
                    }
                    mpz_init_set_str((*arr)[*count], token, 10);
                    (*count)++;
                }
            }
            token = strtok(NULL, ",\n");
        }
    }
    fclose(fp);
}

/* Producto en árbol balanceado in-place: O(N log N) en el tamaño del producto
 * final (vs O(N²) del left-fold). Aprovecha la mul FFT de GMP en niveles altos. */
static void tree_product(mpz_t *arr, size_t count, mpz_t result) {
    if (count == 0) { mpz_set_ui(result, 1); return; }
    while (count > 1) {
        size_t pairs = count / 2;
        for (size_t i = 0; i < pairs; i++) {
            /* Aliasing seguro: dst=arr[i] coincide con src=arr[2i] solo en i=0
             * (mpz_mul soporta dst==src). Para i>0, los índices fuente 2i,2i+1
             * están adelantados respecto a i, así que no se corrompen lecturas. */
            mpz_mul(arr[i], arr[2 * i], arr[2 * i + 1]);
        }
        if (count & 1) {
            mpz_swap(arr[pairs], arr[count - 1]);
            count = pairs + 1;
        } else {
            count = pairs;
        }
    }
    mpz_set(result, arr[0]);
}

int main(int argc, char *argv[]) {

    mpz_t n;
    mpz_init(n);
    mpz_set_str(n, argv[1], 10);
    printf("N:");
    mpz_out_str(stdout, 10, n);
    printf("\n");

    mpz_t *q_arr = NULL;
    size_t q_count = 0, q_cap = 0;

    load_mpz_file("salida.txt", &q_arr, &q_count, &q_cap);
    if (q_count == 0) {
        fprintf(stderr, "Falta archivo salida.txt o está vacío");
        exit(EXIT_FAILURE);
    }
    load_rootas_squared("roota_list.txt", &q_arr, &q_count, &q_cap);

    mpz_t qx;
    mpz_init(qx);
    tree_product(q_arr, q_count, qx);

    for (size_t i = 0; i < q_count; i++) mpz_clear(q_arr[i]);
    free(q_arr);

    if (mpz_sgn(qx) < 0)
        mpz_mul_si(qx, qx, -1);

    /* Se omite mpz_perfect_square_p(qx): el solver GF(2) garantiza que qx es
     * cuadrado. Si fallara, mpz_sqrt daría floor(sqrt) y los GCDs serían
     * triviales — equivalente al "no encontrado" actual. */

    FILE *file = fopen("pos.txt", "r");
    if (file == NULL) {
        fprintf(stderr, "Falta archivo pos.txt");
        exit(EXIT_FAILURE);
    }

    mpz_t mulX, tmp;
    mpz_inits(mulX, tmp, NULL);
    mpz_set_ui(mulX, 1);

    /* mulX = ∏ lhs_i mod N. Solo se usa en gcd(±sqrt(qx) − mulX, N),
     * así que reducimos cada vez y nunca crece más allá de |N|. */
    char buf[BUFSIZ];
    while (fgets(buf, sizeof(buf), file) != NULL) {
        if (buf[0] == '\n' || buf[0] == '\0') continue;
        mpz_set_str(tmp, buf, 10);
        mpz_mul(mulX, mulX, tmp);
        mpz_mod(mulX, mulX, n);
    }
    fclose(file);
    mpz_clear(tmp);

    mpz_t sqr, gcd, res;
    mpz_inits(sqr, gcd, res, NULL);

    mpz_sqrt(sqr, qx);

    mpz_t p, q;
    mpz_inits(p, q, NULL);

    /* Multiplicador Knuth-Schroeppel: factorizamos kN, así que hay que
     * quitar k de los GCDs antes de devolver los factores reales. */
    unsigned long ks_mult = 1;
    {
        FILE *fmult = fopen("multiplier.txt", "r");
        if (fmult) {
            char mbuf[64];
            if (fgets(mbuf, sizeof(mbuf), fmult))
                ks_mult = strtoul(mbuf, NULL, 10);
            fclose(fmult);
        }
    }

    /* gcd(sqr − mulX, N) → p */
    mpz_sub(res, sqr, mulX);
    mpz_gcd(gcd, res, n);
    mpz_set(p, gcd);

    if (ks_mult > 1) {
        mpz_t g_tmp;
        mpz_init(g_tmp);
        mpz_gcd_ui(g_tmp, p, ks_mult);
        while (mpz_cmp_ui(g_tmp, 1) > 0) {
            mpz_divexact(p, p, g_tmp);
            mpz_gcd_ui(g_tmp, p, ks_mult);
        }
        mpz_clear(g_tmp);
        if (mpz_cmp_ui(p, 1) == 0) {
            mpz_clears(p, q, mulX, sqr, gcd, res, qx, n, NULL);
            exit(EXIT_FAILURE);
        }
    }

    /* gcd(sqr + mulX, N) → q */
    mpz_add(res, sqr, mulX);
    mpz_gcd(gcd, res, n);
    mpz_set(q, gcd);

    if (ks_mult > 1) {
        mpz_t g_tmp;
        mpz_init(g_tmp);
        mpz_gcd_ui(g_tmp, q, ks_mult);
        while (mpz_cmp_ui(g_tmp, 1) > 0) {
            mpz_divexact(q, q, g_tmp);
            mpz_gcd_ui(g_tmp, q, ks_mult);
        }
        mpz_clear(g_tmp);
        if (mpz_cmp_ui(q, 1) == 0) {
            mpz_clears(p, q, mulX, sqr, gcd, res, qx, n, NULL);
            exit(EXIT_FAILURE);
        }
    }

    if (mpz_cmp_ui(p, 1) != 0 && mpz_cmp_ui(q, 1) != 0) {
        gmp_printf("P:%Zd \nQ:%Zd\n", p, q);
        mpz_clears(p, q, mulX, sqr, gcd, res, qx, n, NULL);
        exit(EXIT_SUCCESS);
    }

    mpz_clears(p, q, mulX, sqr, gcd, res, qx, n, NULL);
    exit(EXIT_FAILURE);
}
