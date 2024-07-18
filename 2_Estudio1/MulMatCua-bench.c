#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h> // Añadir la biblioteca matemática para usar fabs


#define WARMUP 4
#define REPETICIONES 4

// sequential matrix multiplication
void seq_mult(double *A, double *B, double *C_check, int matrix_size)
{
  int i, j, k;
  double sum;
  for (i = 0; i < matrix_size; i++)
  {
    for (j = 0; j < matrix_size; j++)
    {
      sum = 0.0;
      for (k = 0; k < matrix_size; k++)
      {
        sum += A[i * matrix_size + k] * B[k * matrix_size + j]; // A[i][k] * B[k][j]
      }
      C_check[i * matrix_size + j] = sum; // C_check[i][j]
    }
  }
}

// mm
//  A = A matrix; B = B Matrix; C = Result Matrix ; matrix_size (ldm) = matrix dimms (squared); ldm = matrix leading dimension
void mm(double *A, double *B, double *C, int matrix_size, int num_threads_t, int num_filas_f)
{
  int ldm = matrix_size;
  // Set the number of threads
  omp_set_num_threads(num_threads_t);

  int i, j, k, iam, nprocs;
  double sum;
// Threads are created, all threads have iam and nprocs private variables
#pragma omp parallel private(iam, nprocs)
  {
    nprocs = omp_get_num_threads();
    iam = omp_get_thread_num();

// The matrix multiplication is done in parallel, each thread has its own private variables
//  i: row of the matrix, sum: result of the multiplication for one element of C,
// j: column of the matrix, k: column of A and row of B
#pragma omp for private(i, j, k, sum) schedule(static, num_filas_f)
    for (i = 0; i < matrix_size; i++)
    {
      for (j = 0; j < matrix_size; j++)
      {
        sum = 0.0;
        for (k = 0; k < matrix_size; k++)
        {
          sum += A[i * ldm + k] * B[k * ldm + j]; // A[i][k] * B[k][j]
        }
        C[i * ldm + j] = sum; // C[i][j]
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
int comprobar(double *C, double *C_check, int matrix_size)
{
  int i, j;
  double tolerancia = 0.01; // Definir la tolerancia hasta el segundo decimal
  for (i = 0; i < matrix_size; i++)
  {
    for (j = 0; j < matrix_size; j++)
    {
      // Comprobar si la diferencia absoluta es mayor a la tolerancia
      if (fabs(C[i * matrix_size + j] - C_check[i * matrix_size + j]) > tolerancia)
      {
        return -1;
      }
    }
  }
  return 0;
}

int main(int argc, char *argv[])
{
  // Indicar uso argumentos
  if (argc < 4)
  {
    printf("Uso: %s <dimension_matriz_n> <numero_filas_consecutivas_F> <numero_hilos_t>\n", argv[0]);
    return -1;
  }
  // Comprobar que los argumentos son enteros
  if (atoi(argv[1]) == 0 || atoi(argv[2]) == 0 || atoi(argv[3]) == 0)
  {
    printf("Los argumentos deben ser enteros\n");
    return -1;
  }

  // Convertir argumentos
  int n = atoi(argv[1]);
  int num_filas_f = atoi(argv[2]);
  int num_threads_t = atoi(argv[3]);

  // Variables para medir tiempo
  double start, fin, Gflops, avg_time;

  // Reservar memoria para las matrices linealmente
  double *A = (double *)malloc(sizeof(double) * n * n);
  double *B = (double *)malloc(sizeof(double) * n * n);
  double *C = (double *)malloc(sizeof(double) * n * n);

  // Inicializar matrices A y B aquí
  initializealea(A, n * n);
  initializealea(B, n * n);

  //WARM UP
  for (int i = 0; i < WARMUP; i++)
  {
    mm(A, B, C, n, num_threads_t, num_filas_f);
  }
  
  start=omp_get_wtime();
  //Realizar la multiplicación de matrices REPETICIONES veces
  for (int i = 0; i < REPETICIONES; i++)
  {
    mm(A, B, C, n, num_threads_t, num_filas_f);
  }
  fin=omp_get_wtime();
    
  avg_time = (fin - start) / REPETICIONES;
  if(avg_time==0.) {
    printf("No hay suficiente precision\n");
  }
  else {
    //GFLOPS/s=(2N^3−N^2/avg_matmult_time)/1e9
    Gflops= ((2.*n*n*n-n*n)/(avg_time))/1000000000.;
    //Formato Impresion: <dimension_matriz_n>, <numero_filas_consecutivas_F>, <numero_hilos_t>, <tiempo_medio>, <Gflops>
    printf("%d, %d, %d, %.10lf, %.10lf\n", n, num_filas_f, num_threads_t, avg_time, Gflops);
  }

  // Liberar memoria
  free(A);
  free(B);
  free(C);

  return 0;
}
