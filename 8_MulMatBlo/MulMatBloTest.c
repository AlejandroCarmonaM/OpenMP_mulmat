#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h> // Añadir la biblioteca matemática para usar fabs


// sequential rectangular matrix multiplication
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
// num_hilos = número de hilos
void mm(double *A, double *B, double *C, int m, int n, int k, int tam_blo_b, int num_hilos)
{
    int ldm = n;
    int fila, col, fila_bloque, col_bloque, l;
    double sum;

    int iam;
    omp_set_num_threads(num_hilos);
    //private -> cada hilo tiene su propia copia de la variable
    //collapse(2) -> colapsa los 2 primeros bucles anidados en un solo bucle
    //schedule(static, 1) -> divide el trabajo en bloques de tamaño 1
    #pragma omp parallel for private(fila, col, fila_bloque, col_bloque, l, sum) collapse(2) schedule(static, 1)
    for (fila = 0; fila < m; fila += tam_blo_b)
    {
        for (col = 0; col < k; col += tam_blo_b)
        {
            #if defined (_OPENMP) 
            iam = omp_get_thread_num();
            #endif
            #ifdef DEBUG
            printf("Soy el Hilo %d haciendo el BLOQUE (%d, %d)\n", iam, fila, col);
            #endif
            for (fila_bloque = fila; fila_bloque < fila + tam_blo_b; fila_bloque++)
            {
                for (col_bloque = col; col_bloque < col + tam_blo_b; col_bloque++)
                {
                    sum = 0.0;
                    for (l = 0; l < ldm; l++)
                    {
                        sum += A[fila_bloque * ldm + l] * B[l * k + col_bloque];
                    }
                    C[fila_bloque * k + col_bloque] = sum;
                }
            }
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
    int array_hilos[] = {1, 4, 8, 16};
    int i, m, n, k, tam_bloque;
    double *A, *B, *C, *C_check;

    int is_error = 0;

    for (int i = 0; i < 4; i++)
    {
        int num_hilos_max = array_hilos[i];

        for (int i = 0; i < array_size; i++)
        {
            m = m_array[i];
            n = n_array[i];
            k = k_array[i];
            tam_bloque = tam_bloque_array[i];
            // printf("m = %d, n = %d, k = %d, tam_bloque = %d\n", m, n, k, tam_bloque);
            A = (double *)malloc(m * n * sizeof(double));
            B = (double *)malloc(n * k * sizeof(double));
            C = (double *)malloc(m * k * sizeof(double));
            C_check = (double *)malloc(m * k * sizeof(double));

            initializealea(A, m * n);
            initializealea(B, n * k);

            // Multiplicación secuencial
            seq_mult(A, B, C_check, m, n, k);

            // Multiplicación en paralelo
            mm(A, B, C, m, n, k, tam_bloque, num_hilos_max);

            // Comprobar si el resultado es correcto

            if (comprobar(C, C_check, m, k) == 0)
            {
                //printf("(%d, %d, %d, %d, %d): OK\n", m, n, k, tam_bloque, num_hilos_max);
            }
            else
            {
                //(m, n, k tam_blo_b): ERROR
                printf("(%d, %d, %d, %d, %d): ERROR\n", m, n, k, tam_bloque, num_hilos_max);
                is_error = 1;
            }

            // Liberar memoria
            free(A);
            free(B);
            free(C);
            free(C_check);
        }
    }
    if (is_error == 0)
    {
        printf("Todas las pruebas han pasado\n");
    }
    return 0;
}