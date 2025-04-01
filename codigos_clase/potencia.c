#include <stdio.h>

int power (int base, int n);
int main()
//prueba la funcion power
{
    int i;
    for (i = 0; i < 10; ++i)
        printf("%d %d %d\n", i, power(2,i), power(-3,i));
    return 0;
}

int power(int base, int n)
//power: eleva la base a la potencia n; n >= 0
{
    int i, p;
    p = 1;
    for (i = 1; i <= n; ++i)
        p = p * base;
    return p;
}