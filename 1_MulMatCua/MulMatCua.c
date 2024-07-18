#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h> // Añadir la biblioteca matemática para usar fabs
//Definir DEBUG para imprimir las matrices
#define DEBUG

//Multiplicación secuencial de matrices
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
//  A = Matriz A, B = Matriz B, C = Matriz C, matrix_size (ldm) = tamaño de la matriz,
// num_threads_t = número de hilos, num_filas_f = número de filas consecutivas
void mm(double *A, double *B, double *C, int matrix_size, int num_threads_t, int num_filas_f)
{
  int ldm = matrix_size;
  // Ajustar el número de hilos
  omp_set_num_threads(num_threads_t);

  int i, j, k, iam, nprocs;
  double sum;
// Los hiños se crean, todos los hilos tienen las variables privadas iam y nprocs
#pragma omp parallel private(iam, nprocs)
  {
    nprocs = omp_get_num_threads();
    iam = omp_get_thread_num();

//La multiplicación de matrices se realiza en paralelo, cada hilo tiene sus propias variables privadas:
//  i: fila de la matriz, sum: resultado de la multiplicación para un elemento de C,
//j: columna de la matriz, k: columna de A y fila de B
#pragma omp for private(i, j, k, sum) schedule(static, num_filas_f)
    for (i = 0; i < matrix_size; i++)
    {
      #ifdef DEBUG
      printf("thread %d fila %d \n", iam, i);
      #endif
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
  // Comprobar que los argumentos son enteros mayores a 0
  if (atoi(argv[1]) <= 0 || atoi(argv[2]) <= 0 || atoi(argv[3]) <= 0)
  {
    printf("Los argumentos deben ser enteros mayores a 0\n");
    return -1;
  }

  // Convertir argumentos
  int n = atoi(argv[1]);
  int num_filas_f = atoi(argv[2]);
  int num_threads_t = atoi(argv[3]);

  // Reservar memoria para las matrices linealmente
  double *A = (double *)malloc(sizeof(double) * n * n);
  double *B = (double *)malloc(sizeof(double) * n * n);
  double *C = (double *)malloc(sizeof(double) * n * n);

  // Inicializar matrices A y B aquí
  initializealea(A, n * n);
  initializealea(B, n * n);

  // Multiplicar matrices
  mm(A, B, C, n, num_threads_t, num_filas_f);

  #ifdef DEBUG
  //Comprobar que el resultado es correcto
  double *C_check = (double *)malloc(sizeof(double) * n * n);
  seq_mult(A, B, C_check, n);
  if(comprobar(C, C_check, n) == 0)
  {
    printf("\nMULTIPLICACION: OK\n");
  }
  else
  {
    printf("\nMULTIPLICACION: ERROR\n");
    printf("Matriz C correcta\n");
    escribir(C_check, n, n, n);
  }
  free(C_check);
  
  #endif
  
  // Imprimir matrices
  printf("\nMatriz A\n");
  escribir(A, n, n, n);
  printf("Matriz B\n");
  escribir(B, n, n, n);
  printf("Matriz C Resultado\n");
  escribir(C, n, n, n);
  

  // Liberar memoria
  free(A);
  free(B);
  free(C);



  return 0;
}
