import numpy as np
from scipy.integrate import quad

# Constantes
En = -0.4
Hop_t = -0.2
N = 50  # Asegurarse de que N no sea mayor que el tamaño de A o B
LowerLimitDx = np.pi - np.arccos(En / (4 * Hop_t))
Delta0 = 24.0 * np.sqrt(2.0)

# Leer datos desde archivos
def load_data(file_name):
    A = []
    B = []
    with open(file_name, 'r') as f:
        for line in f:
            # Ignorar líneas vacías
            if not line.strip():
                continue
            try:
                # Intenta descomponer la línea en dos valores
                a, b = map(float, line.split())
                A.append(a)
                B.append(b)
            except ValueError:
                # Si no hay exactamente dos valores en la línea, ignórala
                print(f"Advertencia: línea malformada encontrada en el archivo {file_name}: {line}")
                continue
    return np.array(A), np.array(B)


# Funciones para el cálculo de DOS
def func1(x, n, A, B):
    # Aseguramos que n esté en el rango adecuado
    n = min(n, len(A) - 1)  # Asegura que n no sea mayor que 49 si el tamaño de A es 50
    n = max(n, 0)           # Asegura que n no sea menor que 0
    
    cos_x = np.cos(x)
    sin_x = np.sqrt(1 - cos_x ** 2)
    aux1 = En + 2 * Hop_t * cos_x
    aux2 = Hop_t
    z = np.pi - np.arccos(0.5 * aux1 / aux2)
    cos_z = np.cos(z)
    sin_z = np.sqrt(1 - cos_z ** 2)
    vx = 2 * Hop_t * sin_x
    vz = 2 * Hop_t * sin_z
    vel = np.sqrt(vx ** 2 + vz ** 2)
    der_zx = (vx / vz) ** 2
    factor = np.sqrt(1 + der_zx)
    jacobian = factor / vel
    
    a_k = A[n] ** 2 - B[n] ** 2 - (cos_x - cos_z) ** 2
    b = 2 * A[n] * B[n]
    rho_k = np.sqrt(a_k ** 2 + b ** 2)
    positive_k = 1 + a_k / rho_k
    aux3 = np.sqrt(2.0 * rho_k)
    aux4 = A[n] / aux3
    result = jacobian * aux4 * np.sqrt(positive_k)
    
    return result

def func2(x, n, A, B):
    # Verifica que n no sea mayor que el tamaño de A o B
    n = min(n, len(A) - 1)  # Asegura que n no sea mayor que 49 si el tamaño de A es 50
    n = max(n, 0)  # Asegura que n no sea menor que 0
    
    cos_x = np.cos(x)
    sin_x = np.sqrt(1 - cos_x ** 2)
    aux1 = En + 2 * Hop_t * cos_x
    aux2 = Hop_t
    z = np.pi - np.arccos(0.5 * aux1 / aux2)
    cos_z = np.cos(z)
    sin_z = np.sqrt(1 - cos_z ** 2)
    vx = 2 * Hop_t * sin_x
    vz = 2 * Hop_t * sin_z
    vel = np.sqrt(vx ** 2 + vz ** 2)
    der_zx = (vx / vz) ** 2
    factor = np.sqrt(1 + der_zx)
    jacobian = factor / vel
    
    a_k = A[n] ** 2 - B[n] ** 2 - (cos_x - cos_z) ** 2
    b = 2 * A[n] * B[n]
    rho_k = np.sqrt(a_k ** 2 + b ** 2)
    positive_k = 1 + a_k / rho_k
    aux3 = np.sqrt(2.0 * rho_k)
    aux4 = A[n] / aux3
    result = jacobian * aux4 * np.sqrt(positive_k)
    
    return result

def dos_f(x):
    cos_x = np.cos(x)
    sin_x = np.sqrt(1 - cos_x**2)
    aux1 = En + 2 * Hop_t * cos_x
    z = np.pi - np.arccos(0.5 * aux1 / Hop_t)
    cos_z = np.cos(z)
    sin_z = np.sqrt(1 - cos_z**2)
    vx = 2 * Hop_t * sin_x
    vz = 2 * Hop_t * sin_z
    vel = np.sqrt(vx**2 + vz**2)
    der_zx = (vx / vz)**2
    factor = np.sqrt(1 + der_zx)
    jacobian = factor / vel
    return jacobian

# Función de integración Romberg
def romberg_integration(func, a, b, args=(), tol=1e-6):
    result, _ = quad(func, a, b, args=args)
    return result

# Función principal para procesar los datos
def process_data():
    file_names = ["data/c00g000E04t02.dat", "data/c00g001E04t02.dat", 
                  "data/c00g005E04t02.dat", "data/c00g010E04t02.dat", 
                  "data/c00g015E04t02.dat"]
    
    dos_names = ["output/output_000.dat", "output/output_001.dat", 
                 "output/output_005.dat", "output/output_010.dat", 
                 "output/output_015.dat"]

    for i in range(5):
        A, B = load_data(file_names[i])

        with open(dos_names[i], "w") as dos_file:
            for n in range(min(N, len(A))):  # Asegúrate de que n no supere el tamaño de A o B
                normalization = romberg_integration(dos_f, LowerLimitDx, np.pi)
                f = romberg_integration(func1, LowerLimitDx, np.pi, args=(n, A, B))
                g = romberg_integration(func2, LowerLimitDx, np.pi, args=(n, A, B))
                DOS = (f + g) / normalization
                dos_file.write(f"{A[n]}\t{DOS}\n")

if __name__ == "__main__":
    process_data()
