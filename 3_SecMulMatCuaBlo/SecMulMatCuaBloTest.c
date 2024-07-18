#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

// Multiplicación secuencial de matrices
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
// tam_blo_b = tamaño del bloque
void mm(double *A, double *B, double *C, int matrix_size, int tam_blo_b)
{
  double sum; // Variable para almacenar el resultado de la multiplicación
  int ldm = matrix_size; // Tamaño de la matriz = ldm (matriz cuadrada)
  int fila, col, fila_bloque, col_bloque, k; // Variables para los bucles

  //Avanzamos de tamaño de bloque en tamaño de bloque para recorrer la matriz
  //así en el 3er bucle for se opera el bloque concreto Pos Inicio:(fila, fila+tam_blo_b), Pos fin:(col, col+tam_blo_b)
  for (fila = 0; fila < ldm; fila += tam_blo_b) {
    for (col = 0; col < ldm; col += tam_blo_b) {
      #ifdef DEBUG
      printf("BLOQUE (%d, %d)\n", fila, col); // Imprimir el bloque actual con su Pos Inicio
      #endif
      // Hacemos la multiplicación de matrices para bloques
      for (fila_bloque = fila; fila_bloque < fila+tam_blo_b; fila_bloque ++) {
        for(col_bloque = col; col_bloque < col+tam_blo_b; col_bloque++) {
          sum = 0.0;
          //recorremos toda la fila y la columna para hacer la obtener el resultado de 1 elemento del bloque
          for (k = 0; k < ldm; k++) {
            sum += A[fila_bloque * ldm + k] * B[k * ldm + col_bloque];
          }
          C[fila_bloque * ldm + col_bloque] = sum;
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

    for (int n = n_min; n <= n_max; n++)
    {
        for (int tam_blo_b = 1; tam_blo_b <= n; tam_blo_b++)
        {
            // comprobamos que sea múltiplo
            if (n % tam_blo_b == 0)
            {

                // printf("Procesando matriz de tamaño %d con bloque de tamaño %d\n", n, tam_blo_b);

                // Reserva y inicialización de matrices...
                double *A = (double *)malloc(sizeof(double) * n * n);
                double *B = (double *)malloc(sizeof(double) * n * n);
                double *C = (double *)malloc(sizeof(double) * n * n);
                initializealea(A, n * n);
                initializealea(B, n * n);

                // Ejecución de la multiplicación de matrices...
                mm(A, B, C, n, tam_blo_b);

                // Comprobación de resultados...
                double *C_check = (double *)malloc(sizeof(double) * n * n);
                seq_mult(A, B, C_check, n);
                if (comprobar(C, C_check, n) == 0)
                {
                    //(n, tam_blo_b): OK
                    // printf("(%d, %d):   OK\n", n, tam_blo_b);
                }
                else
                {
                    //(n, tam_blo_b): ERROR
                    printf("(%d, %d):   ERROR\n", n, tam_blo_b);
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
    if (is_error == 0)
    {
        printf("Todas las pruebas han pasado\n");
    }
    return 0;
}
