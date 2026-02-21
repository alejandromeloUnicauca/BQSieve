import sys
import subprocess
import os
import shutil
import time
import struct

PATH_TMP="/tmp/bqsieve"
PATH_BWC="/home/debian/Descargas/cado-nfs/build/debian12/linalg/bwc"
MATRIX_BIN="matrix.bin"
MATRIX_RW_BIN="matrix.rw.bin"
MATRIX_CW_BIN="matrix.cw.bin"

def main():
    remove_temp_files()

    # Detectar -v (verbose) y separarlo de los args que van a B_QSieve
    verbose = "-v" in sys.argv
    bqs_args = sys.argv[1:]  # pasar todos los args incluyendo -v a B_QSieve

    # Parsear el número N de los argumentos
    num = None
    argv = sys.argv[1:]
    i = 0
    while i < len(argv):
        if argv[i] == "-d" and i + 1 < len(argv):
            num = int(argv[i + 1])
            i += 2
        elif argv[i] == "-h" and i + 1 < len(argv):
            num = int(argv[i + 1], 16)
            i += 2
        else:
            i += 1

    if num is None:
        print("Error: se requiere -d <N> o -h <N>", file=sys.stderr)
        sys.exit(1)

    args = ["./B_QSieve"] + bqs_args

    inicioBQS = time.time()
    if verbose:
        exit_status = subprocess.run(args).returncode
    else:
        exit_status = subprocess.run(args, stdout=subprocess.DEVNULL).returncode
    finBQS = time.time()
    tiempoBQS = finBQS - inicioBQS
    if verbose:
        print(f'\nCribado y construccion de matriz:{tiempoBQS}s, exit_status: {exit_status}')
    if exit_status != 0:
        remove_temp_files()
        sys.exit(1)
    
    ascii_to_binary_matrix("matrix.txt",MATRIX_BIN)

    inicioCado = time.time()
    _devnull = subprocess.DEVNULL if not verbose else None
    subprocess.run([PATH_BWC+"/mf_scan2",MATRIX_BIN],
                   stdout=_devnull, stderr=_devnull)

    if os.path.exists(PATH_TMP):
        shutil.rmtree(PATH_TMP)

    os.mkdir(PATH_TMP)
    subprocess.run(["cp", MATRIX_BIN, MATRIX_RW_BIN, MATRIX_CW_BIN, PATH_TMP],
                   stdout=_devnull, stderr=_devnull)

    if verbose:
        print("Solucionando Matriz...")
    subprocess.run([PATH_BWC+"/bwc.pl",
                    ":complete", "thr=2", "m=64", "n=64", "nullspace=left", "interval=100",
                    "matrix="+PATH_TMP+"/"+MATRIX_BIN, "wdir="+PATH_TMP, "interleaving=0"],
                    stdout=open("outputbwc.txt", "w"),
                    stderr=_devnull)

    # CADO-NFS genera K.sols0-64.X.txt solo para el vector que es
    # nullspace válido (X puede ser 0 o 1). Buscar cuál existe.
    ksol_file = None
    for candidate in sorted(os.listdir(PATH_TMP)):
        if candidate.startswith("K.sols") and candidate.endswith(".txt"):
            ksol_file = candidate
            break
    if ksol_file is None:
        print("Error: no se encontró archivo K.sols*.txt en", PATH_TMP)
        remove_temp_files()
        sys.exit(1)
    if verbose:
        print(f"Archivo de solución: {ksol_file}")
    subprocess.run(["cp", PATH_TMP + "/" + ksol_file, "./K.sols.txt"])
    hex_to_binary("K.sols.txt", "vec.txt")
    finCado = time.time()
    tiempoSolM = finCado-inicioCado
    if verbose:
        print(f'\nSolucion Matriz:{tiempoSolM}s')

    inicioMcd = time.time()
    process_polynomial(num, verbose)

    with open("salidap.txt") as f:
        print(f.read())
    finMcd = time.time()

    tiempomcd = finMcd - inicioMcd
    if verbose:
        print(f'\nBusqueda de solucion:{tiempomcd}s')
    tiempof = tiempoBQS+tiempoSolM+tiempomcd
    print(f'Tiempo final:{tiempof}s')

    remove_temp_files()


def process_polynomial(num, verbose=False):
    with open("vec.txt") as vec_file:
        filas = sum(1 for _ in vec_file)

    if verbose:
        print(filas)

    for i in range(1, 64):
        with open("polinomio.txt") as polinomio_file, open("vec.txt") as vec_file:
            for line_polinomio, line_vec in zip(polinomio_file, vec_file):
                f2 = line_vec[i-1]

                if f2 == "1":
                    parts = line_polinomio.strip().split(";")
                    pos = parts[0].strip()
                    f1 = parts[1].strip()
                    # Si hay tercera columna es roota (modo MPQS)
                    rel_roota = parts[2].strip() if len(parts) > 2 else None

                    with open("salida.txt", "a") as salida_file, \
                         open("pos.txt", "a") as pos_file:
                        salida_file.write(f1 + "\n")
                        pos_file.write(pos + "\n")
                    if rel_roota:
                        with open("roota_list.txt", "a") as rl:
                            rl.write(rel_roota + "\n")

        if not os.path.exists("salida.txt"):
            continue

        mulpoli_args = ["./mulPoli", str(num)]
        exit_status = subprocess.run(mulpoli_args, stdout=open("salidap.txt", "a")).returncode
        if exit_status == 0:
            break

        subprocess.run(["rm", "-f", "salida.txt", "pos.txt", "roota_list.txt"])

def hex_to_binary(file_path, output_file):
    """
    Convierte números hexadecimales en un archivo de entrada a binario y guarda los resultados en un archivo de salida.

    Args:
        file_path (str): Ruta del archivo de entrada que contiene números hexadecimales.
        output_file (str): Ruta del archivo de salida donde se guardarán los números binarios convertidos.

    Returns:
        None
    """
    with open(file_path) as f:
        with open(output_file, "w") as out:
            for line in f:
                line = line.strip().upper()
                bin_num = subprocess.run(["bc", "-lq"], input=f"obase=2; ibase=16; {line}\n", capture_output=True, text=True).stdout.strip().replace("\n", "")
                bin_num_padded = bin_num.zfill(64)  # Rellena con ceros a la izquierda si es necesario
                out.write(bin_num_padded + "\n")

def remove_temp_files():
    files_to_remove = [MATRIX_BIN, MATRIX_RW_BIN, MATRIX_CW_BIN, "matrix.txt", "residuos.txt", "outputbwc.txt", "K.sols.txt", "salidap.txt", "vec.txt", "pos.txt", "salida.txt", "polinomio.txt", "roota.txt", "roota_list.txt", "multiplier.txt"]
    [os.remove(file) for file in files_to_remove if os.path.exists(file)]
    if(os.path.exists(PATH_TMP)):
        shutil.rmtree(PATH_TMP)

def ascii_to_binary_matrix(input_file, output_file):
    """
    Convierte una matriz en formato ASCII a un archivo binario.

    Args:
        input_file: Nombre del archivo ASCII de entrada.
        output_file: Nombre del archivo binario de salida.
    """
    # Abre el archivo ASCII para lectura
    with open(input_file, "r") as input_stream:
        it = iter(input_stream)
        
        # Si tienes un encabezado en la primera línea, lee y descártalo.
        # Elimina esta línea si no hay encabezado.
        header_line = next(it)
        
        # Abre el archivo binario para escritura
        with open(output_file, "wb") as data:
            for line in it:
                # Divide la línea y convierte los valores en enteros
                w, *cols = [int(x) for x in line.strip().split()]
                
                # Asegúrate de que w coincide con el número de columnas especificado
                assert w == len(cols), f"Error: {w} no coincide con {len(cols)} en la línea: {line}"
                
                # Escribe el número de columnas (w) en formato binario
                data.write(struct.pack("<I", w))
                
                # Escribe cada índice de columna (cols) en formato binario
                for j in cols:
                    data.write(struct.pack("<I", j))


if __name__ == "__main__":
    main()
