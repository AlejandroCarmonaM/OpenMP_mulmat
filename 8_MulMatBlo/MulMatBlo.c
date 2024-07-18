#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h> // Añadir la biblioteca matemática para usar fabs
// Definir DEBUG para imprimir las matrices
#define DEBUG

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
    // Indicar uso argumentos: MulMatBlo  <dim_mat_m> <dim_mat_n> <dim_mat_k> <tam_blo_b> <num_hil_t>
        if (argc < 6) {
        printf("Uso: %s  <dim_mat_m> <dim_mat_n> <dim_mat_k> <tam_blo_b> <num_hil_t>\n", argv[0]);
        return -1;
    }

    // Comprobar que los argumentos son enteros positivos
    if (atoi(argv[1]) <= 0 || atoi(argv[2]) <= 0 || atoi(argv[3]) <= 0 || atoi(argv[4]) <= 0 || atoi(argv[5]) <= 0)
    {
        printf("Los argumentos deben ser enteros positivos\n");
        return -1;
    }

    // Convertir argumentos
    // m = filas de A
    // n = columnas de A y filas de B
    // k = columnas de B
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int k = atoi(argv[3]);
    int tam_blo_b = atoi(argv[4]);
    int num_hilos = atoi(argv[5]);

    // Reservar memoria para las matrices linealmente
    double *A = (double *)malloc(sizeof(double) * m * n);
    double *B = (double *)malloc(sizeof(double) * n * k);
    double *C = (double *)malloc(sizeof(double) * m * k);

    // Inicializar matrices A y B aquí
    initializealea(A, m * n);
    initializealea(B, n * k);

    // Multiplicar matrices
    mm(A, B, C, m, n, k, tam_blo_b, num_hilos);
    

#ifdef DEBUG
    // Comprobar que el resultado es correcto
    double *C_check = (double *)malloc(sizeof(double) * m * k);
    seq_mult(A, B, C_check, m, n, k);
    if (comprobar(C, C_check, m, k) == 0)
    {
        printf("\nMULTIPLICACION: OK\n");
    }
    else
    {
        printf("\nMULTIPLICACION: ERROR\n");
        printf("Matriz C correcta\n");
        escribir(C_check, m, k, k);
    }
    free(C_check);

#endif

    // Imprimir matrices
    printf("\nMatriz A\n");
    escribir(A, m, n, n);
    printf("Matriz B\n");
    escribir(B, n, k, k);
    printf("Matriz C Resultado\n");
    escribir(C, m, k, k);

    // Liberar memoria
    free(A);
    free(B);
    free(C);

    return 0;
}