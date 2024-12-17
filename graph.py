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
    # Lista de nombres de archivo de salida
    dos_names = ["output/output_000.dat", "output/output_001.dat", 
                 "output/output_005.dat", "output/output_010.dat", 
                 "output/output_015.dat"]
    
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
    plt.title('Densidad de Estados (DoS) vs. Energía', fontsize=18)
    plt.xlabel('Energía (MeV)', fontsize=18)
    plt.ylabel('Densidad de Estados (DoS)', fontsize=18)
    plt.legend()
    plt.grid(True)
    
    # Mostrar el gráfico
    plt.savefig('images/graph.png', dpi=1000) 
    plt.show()

if __name__ == "__main__":
    plot_dos()
