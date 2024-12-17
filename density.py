import matplotlib.pyplot as plt
import numpy as np

# Función para leer datos de los archivos
def read_data(file_name):
    A = []
    DOS = []
    with open(file_name, 'r') as file:
        for line in file:
            if line.strip():  # Ignorar líneas vacías
                a, dos = map(float, line.split())
                A.append(a)
                DOS.append(dos)
    return np.array(A), np.array(DOS)

# Función para graficar los resultados
def plot_dos():
    # Lista de nombres de archivo de entrada (según la imagen)
    dos_names = ["density/density_output_000.dat", "density/density_output_001.dat", 
                 "density/density_output_005.dat", "density/density_output_010.dat", 
                 "density/density_output_015.dat"]
    
    # Crear una figura para el gráfico
    plt.figure(figsize=(10, 6))

    # Colores y etiquetas para los gráficos
    colors = ['b', 'g', 'r', 'c', 'm']
    labels = [r"$\Gamma^+=0.00$", r"$\Gamma^+=0.01$", r"$\Gamma^+=0.05$", r"$\Gamma^+=0.10$", r"$\Gamma^+=0.15$"]
    
    # Leer y graficar los datos de cada archivo
    for i, dos_file in enumerate(dos_names):
        A, DOS = read_data(dos_file)
        plt.plot(A, DOS, label=labels[i], color=colors[i], linestyle='-', marker='o', markersize=0.5)
    
    # Configuración del gráfico
    plt.xlabel(r'$\omega$(meV)', fontsize=18)
    plt.ylabel(r'Im($\tilde{\omega})$(meV)', fontsize=18)
    plt.legend()
    plt.grid(True)
    
    # Mostrar el gráfico
    plt.show()

if __name__ == "__main__":
    plot_dos()
