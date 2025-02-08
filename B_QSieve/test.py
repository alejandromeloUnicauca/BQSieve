import os
import subprocess

keywords = [
    "Intervalo del polinomio:",
    "tiempo de cribado:",
    "Intervalo despues del cribado:",
    "P:",
    "Q:",
    "Tiempo final:"
]

def get_modulus(bits):
    """Obtiene el módulo de la clave RSA desde un archivo PEM."""
    cmd = ["openssl", "rsa", "-noout", "-modulus", "-in", f"../keys/key{bits}.pem"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        modulus = result.stdout.strip().split('=')[-1]
        return modulus
    else:
        return None

def hex_to_dec(hex_value):
    """Convierte un valor hexadecimal a decimal."""
    return int(hex_value, 16)

def extract_lines_grep(filename, keywords):
    """Usa `grep` para extraer líneas con palabras clave de un archivo."""
    grep_pattern = "|".join(keywords)  # Crea un patrón de búsqueda
    result = subprocess.run(["grep", "-E", grep_pattern, filename], capture_output=True, text=True)
    return result.stdout.strip().split("\n")

def run_b_sieve(key, bits, processors):
    """Ejecuta el script B_Sieve.py con la clave y los parámetros adecuados."""
    cmd = ["python3", "B_Sieve.py", "-d", str(key), "-b", "5", "-c", str(processors)]
    output_file = f"results/result{bits}_{processors}.log"
    with open(output_file, "w") as output:
        subprocess.run(cmd, stdout=output)

    selected_lines = extract_lines_grep(f"results/result{bits}_{processors}.log", keywords)
    with open("results/result.log", "a") as log:  # Abrir el archivo dentro del bucle
        log.write(f"{selected_lines}\n")

def main():
    os.makedirs("results", exist_ok=True)
    
    cont = 1
    for bits in range(135, 176, 10):
        print(f"Prueba N°{cont} de {bits} bits\n")
        with open("results/result.log", "a") as log:
            log.write(f"Prueba N°{cont} de {bits} bits \n")
        
        modulus = get_modulus(bits)
        if not modulus:
            print(f"Error obteniendo el módulo para {bits} bits")
            continue
        
        key = hex_to_dec(modulus)
        print(f"N={key}\n")
        print(f"Digitos de N:{len(str(key))}\n")
        with open("results/result.log", "a") as log:
            log.write(f"N={key}\n")
            log.write(f"Digitos de N:{len(str(key))}\n")
        
        # Prueba con 1, 2, 4, 8 y 16 procesadores
        for processors in [1, 2, 4, 8, 16]:
            print(f"UP={processors}\n")
            with open("results/result.log", "a") as log:  # Abrir el archivo dentro del bucle
                log.write(f"UP={processors}\n")
            run_b_sieve(key, bits, processors)
        
        cont += 1

if __name__ == "__main__":
    main()
