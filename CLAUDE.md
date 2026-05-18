# BQSieve — Cuadratic Sieve MPQS/SIQS

Implementación en C del algoritmo **Multiple Polynomial Quadratic Sieve (MPQS/SIQS)** para factorizar números RSA, con un orquestador en Python. Parte del trabajo de grado en la Universidad del Cauca. Toma referencias de msieve v1.46 (J. Papadopoulos) para optimizaciones de criba y manejo de large primes.

## Build

**Requisitos**: `gcc`, `libgmp-dev` (≥6.2), `libmpfr-dev` (≥4.2), `libomp-dev` (OpenMP)

```bash
cd B_QSieve/
make              # produce: B_QSieve, mulPoli, solve_matrix
make debug        # agrega -ggdb3 extra
make clean
```

`primes.txt` debe estar en el directorio de trabajo al ejecutar cualquier binario.

## Ejecución

```bash
# Pipeline completo (recomendado):
python3 B_Sieve.py -d <N_decimal> [-v]
python3 B_Sieve.py -h <N_hex>     [-v]

# Solo fase de criba (debugging):
./B_QSieve -d <N> [-b NBLOCKS] [-c NCORES] [-v]

# Solo fase de álgebra lineal (debugging):
./solve_matrix matrix.txt K.sols.txt

# Solo fase de raíz cuadrada (debugging):
./mulPoli <N_decimal>    # requiere: salida.txt, pos.txt, [roota_list.txt], [multiplier.txt]
```

Flags de `B_QSieve`:
- `-d <N>` / `-h <N>` — N en decimal o hexadecimal
- `-b <n>` — número de bloques para block-division (0 = trial division puro)
- `-c <n>` — número de threads OpenMP
- `-v` — verbose: imprime progreso detallado

## Arquitectura de archivos

```
D:\BQSieve\
├── B_QSieve/
│   ├── structsqs.h      Todas las estructuras de datos del algoritmo QS
│   ├── B_QSieve.c       Entry point C. Parsea N, elige multiplicador
│   │                    Knuth-Schroeppel, genera base de primos, ejecuta
│   │                    el loop SIQS, escribe matrix.txt y polinomio.txt.
│   ├── sieve.c/.h       Criba logarítmica uint8 en bloques de 32KB (msieve-style).
│   │                    Precomputa sqrt(N) mod p y recíprocos para cada primo.
│   │                    sieve_mpqs() → candidatos x donde la acumulación
│   │                    de log(p) supera el umbral.
│   ├── polynomial.c/.h  Generador de polinomios SIQS. Un solo valor de 'a'
│   │                    produce 2^(s-1) valores de 'b' via iteración Gray code.
│   │                    eval_mpqs_Qx() evalúa Q(x) = ax² + 2bx + c.
│   ├── factoring.c/.h   Trial division con recíprocos precomputados.
│   │                    Large Primes: 1LP (partial relation) y 2LP (double LP
│   │                    con Pollard-rho). Escribe polinomio.txt por relación.
│   ├── mulPoli.c        Fase de raíz cuadrada. Lee salida.txt y pos.txt,
│   │                    multiplica Q(x) hasta obtener cuadrado perfecto,
│   │                    calcula GCD(sqrt(∏Q) ± ∏lhs, N) → P, Q.
│   ├── linalg.h/.c      Solver GF(2) nativo. Eliminación Gaussiana con
│   │                    bitpacking de 64 bits sobre [A^T | I]. Encuentra
│   │                    el nullspace izquierdo de la matriz de exponentes.
│   ├── solve_matrix.c   main() del solver standalone: lee matrix.txt,
│   │                    llama linalg_gauss_gf2, escribe K.sols.txt.
│   ├── B_Sieve.py       Orquestador Python: B_QSieve → solve_matrix → mulPoli.
│   ├── test.py          Benchmark con claves RSA de varios tamaños.
│   ├── makefile         Targets: all, debug, clean.
│   └── primes.txt       Lista de primos (uno por línea). Debe estar en CWD.
├── keys/                Claves RSA de prueba (generadas con genKeys.sh)
└── graphify-out/        Grafo de conocimiento generado por /graphify
```

## Flujo end-to-end

```
B_QSieve.c:
  N  ──► multiplicador k (Knuth-Schroeppel)  ──► kN
     ──► generatePrimesBase(): primos p_i con Legendre(kN, p_i) = 1
     ──► sieve_precompute_roots(): sqrt(kN) mod p_i, recíprocos
     ──► loop SIQS:
           generate_mpqs_poly()  ──► Q(x) = ax²+2bx+c  (Gray code b)
           sieve_mpqs()          ──► candidatos x (suma logp > umbral)
           factoringTrial()/factoringBlocks():
             trialDivisionRecip()  ──► vector exponentes mod 2
             insertarNumero()      ──► fila en mat (GF2 XOR)
             large primes          ──► try_combine_partial()
           cuando n_BSuaves == mat.n_rows: terminar
     ──► imprimirMatriz()  ──► matrix.txt
     ──► polinomio.txt (lhs;Qfile;roota por relación)

solve_matrix:
  matrix.txt  ──► Gauss RREF sobre [A^T | I]  ──► 64 null vectors
             ──► K.sols.txt (hex 64-bit por relación)

B_Sieve.py:
  hex_to_binary(K.sols.txt)  ──► vec.txt (64 chars '0'/'1' por línea)
  process_polynomial():
    Para cada bit k (0..62):
      seleccionar relaciones donde vec[j][k] == '1'
      salida.txt ← valores Qfile; pos.txt ← valores lhs
      mulPoli N  ──► si GCD ∉ {1, N}: imprimir P, Q y terminar
```

## Formatos de archivos intermedios

### `matrix.txt` (escrito por B_QSieve, leído por solve_matrix)
```
n_rows n_cols
count col_idx1 col_idx2 ...     ← fila sparse (índices de columnas con bit=1)
...
```
- `n_rows = base.length + 1 + 64` (relaciones: full + combined partials)
- `n_cols = base.length + 1` (col 0 = bit signo de Q(x); cols 1..base.length = primos)
- La diferencia `n_rows - n_cols = 64` garantiza exactamente 64 vectores en el nullspace

### `polinomio.txt` (escrito por factoring.c, leído por process_polynomial)
```
lhs;Qfile;roota               ← relación full
lhs;Qfile;roota1,roota2       ← combined partial (1LP o 2LP)
```
- `lhs = a·x + b`, `Qfile = lhs² - N = a·Q(x)`, `roota` para la fase sqrt

### `K.sols.txt` (escrito por solve_matrix, leído por hex_to_binary)
```
<hex 16 dígitos>    ← una línea por relación (total n_rows líneas)
```
- El bit `63-k` del hex de la línea j = 1 si la relación j pertenece a la solución k
- Bits 1..63 son activos (63 soluciones); bit 0 nunca es leído por process_polynomial

### `vec.txt` (escrito por hex_to_binary, leído por process_polynomial)
- n_rows líneas de exactamente 64 chars '0'/'1', MSB primero
- `vec[j][k]` = `'1'` ↔ relación j está en solución k

## Tabla de parámetros de criba

Tomada de msieve v1.46. Para N entre dos entradas se interpola linealmente.

| bits de N | fb_size  | large_mult | sieve_size   |
|-----------|----------|------------|--------------|
| 64        | 100      | 40         | 64 KB        |
| 128       | 450      | 40         | 64 KB        |
| 183       | 2,000    | 40         | 64 KB        |
| 200       | 3,000    | 50         | 64 KB        |
| 212       | 5,400    | 50         | 192 KB       |
| 233       | 10,000   | 100        | 192 KB       |
| 249       | 27,000   | 100        | 192 KB       |
| 266       | 50,000   | 100        | 192 KB       |
| 283       | 55,000   | 80         | 192 KB       |
| 298       | 60,000   | 80         | 576 KB       |
| 315       | 80,000   | 150        | 576 KB       |
| 332       | 100,000  | 150        | 576 KB       |
| 348       | 140,000  | 150        | 576 KB       |
| 363       | 210,000  | 150        | 832 KB       |
| 379       | 300,000  | 150        | 1,1 MB       |
| 395       | 400,000  | 150        | 1,3 MB       |
| 415       | 500,000  | 150        | 1,6 MB       |
| 440       | 700,000  | 150        | 2,1 MB       |
| 465       | 900,000  | 150        | 3,2 MB       |
| 490       | 1,100,000| 150        | 4,8 MB       |
| 512       | 1,300,000| 150        | 6,4 MB       |

`large_prime_bound = large_mult × p_max`. 2LP se activa automáticamente para N ≥ 283 bits con fb ≥ 800.

## Convenciones del código

- **Aritmética grande**: `mpz_t` (GMP) para todos los enteros de precisión arbitraria.
- **GF(2) XOR**: `insertarNumero(&mat, fila, col, val)` hace XOR: si hay un 1 y se inserta 1, el resultado es 0. Es la operación `+` en GF(2).
- **Columnas de la matriz**: col 0 = bit de signo de Q(x); cols 1..base.length = exponentes de primos (mod 2); cols extra = factores del coeficiente `a` en SIQS.
- **Threads**: `CORES` (global extern) se establece con `-c`; sieve usa `omp_get_thread_num()`.
- **Verbose**: `VERBOSE` (global extern) activado con `-v`; controla todos los `printf` de progreso.
- **Ownership de `exponents`**: `partial_entry.exponents` es un array allocado; quien acepta la entrada es responsable del `free()`.
- **Multiplier.txt**: si existe, contiene el multiplicador k de Knuth-Schroeppel; `mulPoli` lo lee para eliminar sus factores de p y q antes de imprimir el resultado.

## Solver GF(2): linalg.c

El solver implementa **Montgomery Block Lanczos sobre GF(2)** con tamaño de bloque N = 64 y paralelismo OpenMP:

1. Trabaja sobre la matriz simétrica `B = A × A^T`. El nullspace de B coincide con el nullspace izquierdo de A cuando `rank(A) = n_cols` (caso genérico en MPQS).
2. Mantiene un bloque `V_k` de `n_rows × 64` bits (un `uint64_t` por fila).
3. **Recurrencia de tres términos**:
   - `W_k = B × V_k` (dos multiplicaciones dispersas: `A^T V_k` luego `A W`).
   - `S_k = V_k^T × W_k` (matriz de Gram 64×64).
   - `D_k = S_k^{-1} × (W_k^T W_k)`, `E_k = S_{k-1}^{-1} × (W_{k-1}^T W_k)`.
   - `V_{k+1} = W_k + V_k D_k + V_{k-1} E_k`.
4. En cada iteración se comprueban qué columnas de `A^T V_k` son cero → columna correspondiente de `V_k` es un null vector de `A^T`.
5. Termina cuando se encuentran 64 null vectors o se superan `n_cols/64 × 2 + 256` iteraciones.

Las multiplicaciones dispersas (`matvec_AT`, `matvec_A`) se paralelizan con `#pragma omp parallel for`; las operaciones 64×64 (inner products, rightmul) también usan OpenMP en el bucle sobre filas.

**Complejidad**: O(n_cols/64 × nnz) tiempo; O(n_rows × 8) bytes de memoria.
- n=10,000 (≈200 bits): ~0.01 s, ~1 MB
- n=50,000 (≈300 bits): ~0.3 s, ~6 MB
- n=140,000 (≈348 bits): ~2 s, ~18 MB

## Extensiones futuras

- **Block Lanczos GF(2)** (msieve-style): reemplaza la eliminación gaussiana para N > 350 bits (fb > 140,000). Complejidad O(n × nnz / 64); memoria O(n). Permite escalar hasta 512 bits sin los ~2.5 GB de RAM.
- **Integración directa del solver en B_QSieve.c**: en lugar de escribir `matrix.txt` y releer, llamar `linalg_gauss_gf2()` directamente desde `main()` después de `imprimirMatriz()`, eliminando I/O de disco.
- **Port Windows nativo (MSVC)**: `mulPoli.c` y `sieve.c` usan `<unistd.h>`; reemplazar con guardas `#ifdef _WIN32`. `makefile` necesita `ifeq ($(OS),Windows_NT)` para `del /f` y sufijo `.exe`.
- **Paralelización del loop SIQS**: cada polinomio es independiente; el array de criba por polinomio puede procesarse en threads separados con un mutex sobre la escritura a `mat` y `polinomio.txt`.
