#include <iostream>
#include <fstream>
using namespace std;

/**
 * Lee una matriz desde un archivo
 * @param filename Nombre del archivo
 * @param rows Número de filas
 * @param cols Número de columnas
 * @return Puntero a la matriz creada
 */
float** leerMatrizArchivo(const char* filename, int& rows, int& cols) {
    ifstream archivo(filename);
    if (!archivo.is_open()) {
        throw runtime_error("No se pudo abrir el archivo");
    }

    archivo >> rows >> cols;
    
    float** matriz = new float*[rows];
    for(int i = 0; i < rows; i++) {
        matriz[i] = new float[cols];
    }

    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            archivo >> matriz[i][j];
        }
    }

    archivo.close();
    return matriz;
}

/**
 * Multiplica todos los elementos de la matriz por 2
 * @param matriz Puntero a la matriz
 * @param rows Número de filas
 * @param cols Número de columnas
 */
void multiplicarMatriz(float** matriz, int rows, int cols) {
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            matriz[i][j] *= 2;
        }
    }
}

/**
 * Libera la memoria asignada a la matriz
 * @param matriz Puntero a la matriz
 * @param rows Número de filas
 */
void liberarMatriz(float** matriz, int rows) {
    for(int i = 0; i < rows; i++) {
        delete[] matriz[i];
    }
    delete[] matriz;
}

int main() {
    try {
        int filas, columnas;
        float** matriz;

        // Leer matriz del archivo
        matriz = leerMatrizArchivo("matriz.txt", filas, columnas);

        // Multiplicar elementos por 2
        multiplicarMatriz(matriz, filas, columnas);

        // Mostrar resultado
        cout << "Matriz después de multiplicar por 2:\n";
        for(int i = 0; i < filas; i++) {
            for(int j = 0; j < columnas; j++) {
                cout << matriz[i][j] << " ";
            }
            cout << endl;
        }

        // Liberar memoria
        liberarMatriz(matriz, filas);

    } catch(const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}