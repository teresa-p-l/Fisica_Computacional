#include <iostream>
#include <vector>
#include <limits> // Para std::numeric_limits

const double G = 6.67430e-11; // Constante de gravitación universal en m^3 kg^-1 s^-2

void generateAndMultiplyMatrix(int n, int m) {
    // Crear una matriz n x m con todos los elementos iguales a 1
    std::vector<std::vector<double>> matrix(n, std::vector<double>(m, 1.0));

    // Multiplicar cada elemento por la constante G
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            matrix[i][j] *= G;
        }
    }

    // Imprimir la matriz resultante
    std::cout << "Matriz resultante (" << n << "x" << m << ") multiplicada por G:\n";
    for (const auto& row : matrix) {
        for (const auto& element : row) {
            std::cout << element << " ";
        }
        std::cout << "\n";
    }
}

int main() {
    int n, m;

    // Solicitar al usuario las dimensiones de la matriz con validación
    std::cout << "Introduce el número de filas (n): ";
    std::cin >> n;

    std::cout << "Introduce el número de columnas (m): ";
    std::cin >> m;
    // Generar y multiplicar la matriz
    generateAndMultiplyMatrix(n, m);

    return 0;
}