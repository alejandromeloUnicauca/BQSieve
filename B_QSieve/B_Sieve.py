import sys
import subprocess
import os
import shutil

PATH_TMP="/tmp/bqsieve"
PATH_BWC="/home/debian/Descargas/cado-nfs/build/debian12/linalg/bwc"
MATRIX_BIN="matrix.bin"
MATRIX_RW_BIN="matrix.rw.bin"
MATRIX_CW_BIN="matrix.cw.bin"


def main():
    if len(sys.argv) < 2 or sys.argv[1] == "-h":
        num = int(sys.argv[2], 16)
    else:
        num = int(sys.argv[2])

    exit_status = subprocess.run(["./B_QSieve", sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]]).returncode
    print("exit_status",exit_status)
    if exit_status != 0:
        remove_temp_files()
        sys.exit(1)

    subprocess.run([PATH_BWC+"/mf_scan",
                    "--ascii-in", "mfile=matrix.txt", "--binary-out", "ofile=matrix.bin", "--freq"])

    if os.path.exists(PATH_TMP):
        shutil.rmtree(PATH_TMP)

    os.mkdir(PATH_TMP)
    subprocess.run(["cp", MATRIX_BIN, MATRIX_RW_BIN, MATRIX_CW_BIN, PATH_TMP])

    print("Solucionando Matriz...")
    subprocess.run([PATH_BWC+"/bwc.pl",
                    ":complete", "thr=2", "m=64", "n=64", "nullspace=left", "interval=100",
                    "matrix="+PATH_TMP+"/"+MATRIX_BIN, "wdir="+PATH_TMP, "interleaving=0"],
                    stdout=open("outputbwc.txt", "w"))

    subprocess.run(["cp", PATH_TMP+"/K.sols0-64.0.txt", "/home/debian/prsa/B_QSieve/"])
    hex_to_binary("K.sols0-64.0.txt", "vec.txt")
    process_polynomial(num)

    with open("salidap.txt") as f:
        print(f.read())

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

if __name__ == "__main__":
    main()
