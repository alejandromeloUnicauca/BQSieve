import sys
import subprocess
import os
import time


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

    inicioCado = time.time()
    if verbose:
        print("Solucionando Matriz...")
    result = subprocess.run(
        ["./solve_matrix", "matrix.txt", "K.sols.txt"],
        stdout=subprocess.DEVNULL if not verbose else None,
        stderr=subprocess.DEVNULL if not verbose else None,
    )
    if result.returncode != 0:
        print("Error: solve_matrix falló", file=sys.stderr)
        remove_temp_files()
        sys.exit(1)
    hex_to_binary("K.sols.txt", "vec.txt")
    finCado = time.time()
    tiempoSolM = finCado - inicioCado
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
    tiempof = tiempoBQS + tiempoSolM + tiempomcd
    print(f'Tiempo final:{tiempof}s')

    remove_temp_files()


def process_polynomial(num, verbose=False):
    with open("polinomio.txt") as f:
        relations = []
        for line in f:
            parts = line.strip().split(";")
            if len(parts) < 2:
                continue
            lhs = parts[0].strip()
            qfile = parts[1].strip()
            roota = parts[2].strip() if len(parts) > 2 else None
            relations.append((lhs, qfile, roota))

    with open("vec.txt") as f:
        vecs = [line.rstrip("\n") for line in f]

    if verbose:
        print(len(vecs))

    n = min(len(relations), len(vecs))

    for i in range(1, 64):
        selected = [j for j in range(n) if len(vecs[j]) > i - 1 and vecs[j][i - 1] == "1"]
        if not selected:
            continue

        with open("salida.txt", "w") as fsal, open("pos.txt", "w") as fpos:
            fsal.write("\n".join(relations[j][1] for j in selected) + "\n")
            fpos.write("\n".join(relations[j][0] for j in selected) + "\n")

        rootas = [relations[j][2] for j in selected if relations[j][2] is not None]
        if rootas:
            with open("roota_list.txt", "w") as frl:
                frl.write("\n".join(rootas) + "\n")

        mulpoli_args = ["./mulPoli", str(num)]
        exit_status = subprocess.run(
            mulpoli_args, stdout=open("salidap.txt", "a")
        ).returncode
        if exit_status == 0:
            break

        for f in ["salida.txt", "pos.txt", "roota_list.txt"]:
            if os.path.exists(f):
                os.remove(f)


def hex_to_binary(file_path, output_file):
    """Convierte hex de 16 dígitos a cadenas binarias de 64 chars (MSB primero)."""
    with open(file_path) as f, open(output_file, "w") as out:
        for line in f:
            line = line.strip()
            if not line:
                continue
            val = int(line, 16)
            out.write(format(val, '064b') + "\n")


def remove_temp_files():
    files_to_remove = [
        "matrix.txt", "K.sols.txt", "vec.txt",
        "salidap.txt", "pos.txt", "salida.txt",
        "polinomio.txt", "roota.txt", "roota_list.txt",
        "multiplier.txt", "residuos.txt",
    ]
    for f in files_to_remove:
        if os.path.exists(f):
            os.remove(f)


if __name__ == "__main__":
    main()
