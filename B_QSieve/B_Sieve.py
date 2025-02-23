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
    if len(sys.argv) > 2 and sys.argv[1] == "-h":
        num = int(sys.argv[2], 16)
    elif len(sys.argv) > 2:
        num = int(sys.argv[2])

    args = ["./B_QSieve"] + sys.argv[1:]

    inicioBQS = time.time()
    exit_status = subprocess.run(args).returncode
    finBQS = time.time()
    tiempoBQS = finBQS - inicioBQS
    print(f'\nCribado y construccion de matriz:{tiempoBQS}s, exit_status: {exit_status}')
    if exit_status != 0:
        remove_temp_files()
        sys.exit(1)
    
    ascii_to_binary_matrix("matrix.txt",MATRIX_BIN)

    inicioCado = time.time()
    subprocess.run([PATH_BWC+"/mf_scan2",MATRIX_BIN])

    if os.path.exists(PATH_TMP):
        shutil.rmtree(PATH_TMP)

    os.mkdir(PATH_TMP)
    subprocess.run(["cp", MATRIX_BIN, MATRIX_RW_BIN, MATRIX_CW_BIN, PATH_TMP])

    print("Solucionando Matriz...")
    subprocess.run([PATH_BWC+"/bwc.pl",
                    ":complete", "thr=2", "m=64", "n=64", "nullspace=left", "interval=100",
                    "matrix="+PATH_TMP+"/"+MATRIX_BIN, "wdir="+PATH_TMP, "interleaving=0"],
                    stdout=open("outputbwc.txt", "w"))

    subprocess.run(["cp", PATH_TMP+"/K.sols0-64.0.txt", "./"])
    hex_to_binary("K.sols0-64.0.txt", "vec.txt")
    finCado = time.time()
    tiempoSolM = finCado-inicioCado
    print(f'\nSolucion Matriz:{tiempoSolM}s')

    inicioMcd = time.time()
    process_polynomial(num)

    with open("salidap.txt") as f:
        print(f.read())
    finMcd = time.time()

    tiempomcd = finMcd - inicioMcd
    print(f'\nBusqueda de solucion:{tiempomcd}s')
    tiempof = tiempoBQS+tiempoSolM+tiempomcd
    print(f'Tiempo final:{tiempof}s')

    remove_temp_files()


def process_polynomial(num):
    with open("vec.txt") as vec_file:
        filas = sum(1 for _ in vec_file)

    print(filas)

    for i in range(1, 64):
        with open("polinomio.txt") as polinomio_file, open("vec.txt") as vec_file:
            for line_polinomio, line_vec in zip(polinomio_file, vec_file):
                f2 = line_vec[i-1]

                if f2 == "1":
                    f1 = line_polinomio.split(";")[1].strip()
                    pos = line_polinomio.split(";")[0].strip()

                    with open("salida.txt", "a") as salida_file, open("pos.txt", "a") as pos_file:
                        salida_file.write(f1 + "\n")
                        pos_file.write(pos + "\n")

        exit_status = subprocess.run(["./mulPoli", str(num)], stdout=open("salidap.txt", "a")).returncode
        if exit_status == 0:
            break

        subprocess.run(["rm", "salida.txt", "pos.txt"])

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
    files_to_remove = [MATRIX_BIN, MATRIX_RW_BIN, MATRIX_CW_BIN, "matrix.txt", "residuos.txt", "outputbwc.txt", "K.sols0-64.0.txt", "salidap.txt", "vec.txt", "pos.txt", "salida.txt", "polinomio.txt"]
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
