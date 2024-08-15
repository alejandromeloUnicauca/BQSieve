# B_Qsieve

B_QSieve es la implementación del algoritmo Quadratic sieve mejorado al hacer divisiones por bloques en lugar de divisiones triviales, lo que permite una aceleración sobre el paso del cribado.

## Installation

1. Descargar e instale la librería [GMP 6.2.1](https://gmplib.org/list-archives/gmp-announce/2020-November/000049.html).
```console
./configure
make
make check
make install
```
2. Descargar e instalar la librería [MPFR 4.2.1](https://www.mpfr.org/mpfr-current/#download).

```console
./configure
make
make check
make install
```
3. Instalar [Cado NFS](https://gitlab.inria.fr/cado-nfs/cado-nfs) para la solución de matrices
4. Configurar PATH´s en B_Sieve.py
    - configurar PATH_TMP con una ruta temporal para almacenar las matrices. Ejm: /tmp/bqsieve/
    - Configurar PATH_BWC con la ruta de bwc donde se instalo CADO-NFS. Ejm: /CADO-NFS/build/SO/linalg/bwc
5. Compilar el algoritmo
```console
Make
```

## Usage
```console
Python B_Sieve.py (-d | -h) <N> [-b <NBLOCKS>]
Opciones:
-d # Especifica que el numero N es decimal
-h # Especifica que el numero N es hexadecimal
-b # Al usar esta opcion se deben especificar el numero de Bloques
```
## License
