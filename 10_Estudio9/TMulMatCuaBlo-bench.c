#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h> // Añadir la biblioteca matemática para usar fabs
// Definir DEBUG para imprimir las matrices
#define DEBUG


#define WARMUP 4
#define REPETICIONES 4

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

// mm rectangular
//  A = A matrix; B = B Matrix; C = Result Matrix ; 
// mm rectangular
//  A = A matrix; B = B Matrix; C = Result Matrix ; 
void mm(double *A, double *B, double *C, int m, int n, int k, int tam_blo_b, int num_hilos)
{
    int ldm = n;
    int fila, col, fila_bloque, col_bloque, l;
    double sum;

    int iam;
    omp_set_num_threads(num_hilos);

    //TASK PARALLELISM
    
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            for (fila = 0; fila < m; fila += tam_blo_b)
            {
                for (col = 0; col < k; col += tam_blo_b)
                {
                    #pragma omp task firstprivate(fila, col, fila_bloque, col_bloque, l, sum)
                    {
                        iam = omp_get_thread_num();
                        //printf("Soy el Hilo %d haciendo el BLOQUE (%d, %d)\n", iam, fila, col);
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
    // Indicar uso argumentos: SecMulMatCuaBlo <dimension_matriz_n> <tam_blo_b> <num_hil_t>
        if (argc < 4) {
        printf("Uso: %s <dimension_matriz_n> <tam_blo_b> <num_hil_t>\n", argv[0]);
        return -1;
    }

    
    // Variables para medir tiempo
    double start, fin, Gflops, avg_time;

    // Convertir argumentos
    // m = filas de A
    // n = columnas de A y filas de B
    // k = columnas de B
    int m = atoi(argv[1]);
    int n = atoi(argv[1]);
    int k = atoi(argv[1]);
    int tam_blo_b = atoi(argv[2]);
    int num_hilos = atoi(argv[3]);

    // Reservar memoria para las matrices linealmente
    double *A = (double *)malloc(sizeof(double) * m * n);
    double *B = (double *)malloc(sizeof(double) * n * k);
    double *C = (double *)malloc(sizeof(double) * m * k);

    // Inicializar matrices A y B aquí
    initializealea(A, m * n);
    initializealea(B, n * k);

    // Multiplicar matrices
    mm(A, B, C, m, n, k, tam_blo_b, num_hilos);
    // WARM UP
    for (int i = 0; i < WARMUP; i++)
    {
        mm(A, B, C, m, n, k, tam_blo_b, num_hilos);
    }

    start = omp_get_wtime();
    // Realizar la multiplicación de matrices REPETICIONES veces
    for (int i = 0; i < REPETICIONES; i++)
    {
        mm(A, B, C, m, n, k, tam_blo_b, num_hilos);
    }
    fin = omp_get_wtime();

    avg_time = (fin - start) / REPETICIONES;
    if (avg_time == 0.)
    {
        printf("No hay suficiente precision\n");
    }
    else
    {
        // GFLOPS/s=(2N^3−N^2/avg_matmult_time)/1e9
        Gflops = ((2. * n * n * n - n * n) / (avg_time)) / 1000000000.;
        // Formato Impresion: <dimension_matriz_n>, <tamano_bloque_b>, <numero_hilos_t>, <tiempo_medio>, <Gflops>
        printf("%d, %d, %d, %.10lf, %.10lf\n", n, tam_blo_b, num_hilos, avg_time, Gflops);
    }

    // Liberar memoria
    free(A);
    free(B);
    free(C);

    return 0;
}
