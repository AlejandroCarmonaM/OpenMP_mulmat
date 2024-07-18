#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h> // Añadir la biblioteca matemática para usar fabs


//Multiplicación secuencial de matrices
void seq_mult(double *A, double *B, double *C_check, int m, int n, int k)
{
    int i, j, l;
    double sum;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            sum = 0.0;
            for (l = 0; l < n; l++)
            {
                sum += A[i * n + l] * B[l * k + j]; // A[i][l] * B[l][j]
            }
            C_check[i * k + j] = sum; // C_check[i][j]
        }
    }
}

// mm
// A = Matriz A, B = Matriz B, C = Matriz C, m = filas de A, 
// n = columnas de A y filas de B, k = columnas de B y tam_blo_b = tamaño del bloque
void mm(double *A, double *B, double *C, int m, int n, int k, int tam_blo_b)
{
    int ldm = n; // denominamos ldm (leading dimension) 
    //al número de columnas de A y filas de B puesto que será el valor que ussemos para
    //recorrer las filas de A y las columnas de B en su totalidad
    int fila, col, fila_bloque, col_bloque, l; // Variables para los bucles
    double sum; // Variable para almacenar el resultado de la multiplicación

    //Avanzamos de tamaño de bloque en tamaño de bloque para recorrer la matriz
    //así en el 3er bucle for se opera el bloque concreto Pos Inicio:(fila, fila+tam_blo_b), Pos fin:(col, col+tam_blo_b)
    for (fila = 0; fila < m; fila += tam_blo_b) //Avanzamos de bloque en bloque hasta m
    {
        for (col = 0; col < k; col += tam_blo_b) //Avanzamos de bloque en bloque hasta k
        {
            #ifdef DEBUG
            printf("BLOQUE (%d, %d)\n", fila, col); // Imprimir el bloque actual con su Pos Inicio
            #endif
            //Estos dos bucles se mantienen igual que en la multiplicación de matrices cuadradas
            // Hacemos la multiplicación de matrices para bloques
            for (fila_bloque = fila; fila_bloque < fila + tam_blo_b; fila_bloque++)
            {
                for (col_bloque = col; col_bloque < col + tam_blo_b; col_bloque++)
                {
                    sum = 0.0;
                    //recorremos toda la fila y la columna para hacer la obtener el resultado de 1 elemento del bloque
                    for (l = 0; l < ldm; l++)
                    {
                        sum += A[fila_bloque * ldm + l] * B[l * k + col_bloque];
                    }
                    C[fila_bloque * k + col_bloque] = sum;
                    #ifdef DEBUG
                    printf("%.4lf ", sum); //Imprimimos el valor de un elemento del bloque
                    #endif
                }
                #ifdef DEBUG
                printf("\n"); //Salto de línea para separar las filas del bloque
                #endif
            }
            #ifdef DEBUG
            printf("\n"); //Salto de línea para separar los bloques
            #endif
        }

    }

}


///////////////////////////////////////////////////////////////

void initialize(double *m, int t)
{
    int i;
    for (i = 0; i < t; i++)
        m[i] = (double)(i);
}

///////////////////////////////////////////////////////////////
void initializealea(double *m, int t)
{
    int i;
    for (i = 0; i < t; i++)
        m[i] = (double)rand() / RAND_MAX;
}
///////////////////////////////////////////////////////////////
void escribir(double *m, int fm, int cm, int ldm)
{
    int i, j;
    for (i = 0; i < fm; i++)
    {
        for (j = 0; j < cm; j++)
            printf("%.4lf ", m[i * ldm + j]);
        printf("\n");
    }
}

///////////////////////////////////////////////////////////////
int comprobar(double *C, double *C_check, int m, int k)
{
    int i, j;
    double tolerancia = 0.01; // Definir la tolerancia hasta el segundo decimal
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            // Comprobar si la diferencia absoluta es mayor a la tolerancia
            if (fabs(C[i * k + j] - C_check[i * k + j]) > tolerancia)
            {
                return -1;
            }
        }
    }
    return 0;
}

int main(int argc, char *argv[])
{
    int array_size = 7;
    int m_array[] = {8, 16, 16, 16, 24, 24, 25};
    int n_array[] = {12, 24, 24, 24, 32, 32, 50};
    int k_array[] = {16, 32, 32, 32, 40, 40, 15};
    int tam_bloque_array[] = {4, 4, 8, 2, 4, 8, 5};
    int i, m, n, k, tam_bloque;
    double *A, *B, *C, *C_check;

    int is_error=0;
    for(int i=0; i< array_size; i++)
    {
        m = m_array[i];
        n = n_array[i];
        k = k_array[i];
        tam_bloque = tam_bloque_array[i];
        //printf("m = %d, n = %d, k = %d, tam_bloque = %d\n", m, n, k, tam_bloque);
        A = (double *)malloc(m * n * sizeof(double));
        B = (double *)malloc(n * k * sizeof(double));
        C = (double *)malloc(m * k * sizeof(double));
        C_check = (double *)malloc(m * k * sizeof(double));

        initializealea(A, m * n);
        initializealea(B, n * k);

        // Multiplicación secuencial
        seq_mult(A, B, C_check, m, n, k);

        // Multiplicación en paralelo
        mm(A, B, C, m, n, k, tam_bloque);

        // Comprobar si el resultado es correcto

        if (comprobar(C, C_check, m, k) == 0)
        {
            printf("(%d, %d, %d, %d): OK\n", m, n, k, tam_bloque);
        }
        else
        {
            //(m, n, k tam_blo_b): ERROR
            printf("(%d, %d, %d, %d): ERROR\n", m, n, k, tam_bloque);
            is_error = 1;
        }

        // Liberar memoria
        free(A);
        free(B);
        free(C);
        free(C_check);
    }
    if (is_error == 0)
    {
        printf("Todas las pruebas han pasado\n");
    }
    return 0;
}
