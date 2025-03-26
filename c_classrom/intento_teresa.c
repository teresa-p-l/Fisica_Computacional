#include <stdio.h>

int main() {
    int num1, num2, sum;

    // Solicitar al usuario que ingrese dos números
    printf("Ingrese el primer número: ");
    scanf("%d", &num1);

    printf("Ingrese el segundo número: ");
    scanf("%d", &num2);

    // Calcular la suma de los dos números
    sum = num1 + num2;

    // Mostrar el resultado
    printf("La suma de %d y %d es: %d\n", num1, num2, sum);

    return 0;
}