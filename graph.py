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
    labels = ["g=0.0", "g=0.1", "g=0.5", "g=1.0", "g=1.5"]
    
    # Leer y graficar los datos de cada archivo
    for i, dos_file in enumerate(dos_names):
        A, DOS = read_data(dos_file)
        plt.plot(A, DOS, label=labels[i], color=colors[i], linestyle='-', marker='o')
    
    # Configuración del gráfico
    plt.title('Densidad de Estados (DOS) vs. Energía')
    plt.xlabel('Energía (MeV)')
    plt.ylabel('Densidad de Estados (DOS)')
    plt.legend()
    plt.grid(True)
    
    # Mostrar el gráfico
    plt.show()

if __name__ == "__main__":
    plot_dos()
