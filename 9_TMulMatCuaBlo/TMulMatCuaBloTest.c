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

// A = Matriz A, B = Matriz B, C = Matriz C, m = filas de A, 
// n = columnas de A y filas de B, k = columnas de B y tam_blo_b = tamaño del bloque 
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
        /*single indica que el bloque de código siguiente debe ser ejecutado por un solo hilo.
        nowait sirve para que los otros hilos no necesiten esperar a que el hilo de single termine,
        pueden continuar con la ejecución del código que sigue después del bloque single.*/
        #pragma omp single nowait
        {
            for (fila = 0; fila < m; fila += tam_blo_b)
            {
                for (col = 0; col < k; col += tam_blo_b)
                {
                    /*task firstprivate define una tarea paralela (task) donde las variables
                    fila, col, fila_bloque, col_bloque, l, sum son privadas (firstprivate) para cada tarea.
                    firstprivate garantiza que el valor de las variables para cada hilo será*/
                    #pragma omp task firstprivate(fila, col, fila_bloque, col_bloque, l, sum)
                    {
                        #if defined (_OPENMP) 
                        iam = omp_get_thread_num();
                        #endif
                        #ifdef DEBUG
                        printf("Tarea creada para Hilo %d y BLOQUE (%d, %d)\n", iam, fila, col);
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
    if (argc < 3)
    {
        printf("Uso: %s <n_min> <n_max>\n", argv[0]);
        return -1;
    }

    int n_min = atoi(argv[1]);
    int n_max = atoi(argv[2]);
    int is_error = 0;

    if (n_min <= 0 || n_max < n_min)
    {
        printf("Argumentos no válidos.\n");
        return -1;
    }

    for(int num_hilos = 1; num_hilos <= 16; num_hilos++){
        for (int n = n_min; n <= n_max; n++)
        {
            for (int tam_blo_b = 1; tam_blo_b <= n; tam_blo_b++)
            {
                // comprobamos que sea múltiplo
                if (n % tam_blo_b == 0)
                {

                    //printf("Procesando matriz de tamaño %d con bloque de tamaño %d\n", n, tam_blo_b);

                    // Reserva y inicialización de matrices...
                    double *A = (double *)malloc(sizeof(double) * n * n);
                    double *B = (double *)malloc(sizeof(double) * n * n);
                    double *C = (double *)malloc(sizeof(double) * n * n);
                    initializealea(A, n * n);
                    initializealea(B, n * n);

                    // Ejecución de la multiplicación de matrices...
                    mm(A, B, C, n, n, n, tam_blo_b, num_hilos);

                    // Comprobación de resultados...
                    double *C_check = (double *)malloc(sizeof(double) * n * n);
                    seq_mult(A, B, C_check, n, n, n);
                    if (comprobar(C, C_check, n, n) == 0)
                    {
                        //(n, tam_blo_b): OK
                        //printf("(%d, %d):   OK\n", n, tam_blo_b, num_hilos);
                    }
                    else
                    {
                        //(n, tam_blo_b): ERROR
                        printf("(%d, %d, %d):   ERROR\n", n, tam_blo_b, num_hilos);
                        is_error = 1;
                        
                    }

                    // Limpieza...
                    free(A);
                    free(B);
                    free(C);
                    free(C_check);
                }
            }
        }
    }
    if(is_error == 0) {
        printf("Todas las pruebas han pasado\n");
    }
    return 0;
}
