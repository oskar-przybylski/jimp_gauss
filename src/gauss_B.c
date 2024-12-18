#include "gauss_B.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#define DEBUG 0
/**
 * Zwraca 0 - elimnacja zakonczona sukcesem
 * Zwraca 1 - macierz osobliwa - dzielenie przez 0
 */
void mat_assert(Matrix *mat, const char *func_name) {
  if (DEBUG)
    fprintf(stderr, "\nFunction that called %s: %s \n\n", __func__, func_name);
  assert(mat != NULL && mat->data != NULL &&
         "Podany wskaznik na macierz jest pusty(NULL) lub mat->data jest "
         "puste(NULL)! (mat_assert)");
  assert(mat->r != 0 && mat->c != 0);
  for (int i = 0; i < mat->r; i++) {
    assert(mat->data[i] != NULL &&
           "W macierzy jedna lub wiecej kolumn jest pusta(NULL)! (mat_assert)");
  }
}

void mat_assert_square(Matrix *mat, const char *func_name) {
  if (DEBUG)
    fprintf(stderr, "\nFunction that called %s: %s \n\n", __func__, func_name);
  mat_assert(mat, __func__);
  assert(mat->r == mat->c &&
         "Macierz nie jest kwadratowa! (mat_assert_square)");
}

void dot_by_scalar(
    Matrix *mat, int r,
    double scalar) { // indeksuje od 0; r = 0 to pierwszy wiersz macierzy
  mat_assert(mat, __func__);
  for (int i = 0; i < mat->c; i++) {
    mat->data[r][i] *= scalar;
  }
}

void rows_substract(
    Matrix *mat, int r1,
    int r2) { // r1 to wiersz modyfikowany. od wiersza r1 odejmujemy r2
  mat_assert(mat, __func__);
  assert(r1 < mat->r && "wiersz 1 (r1) jest jest wiekszy niz ilosc wierszow w "
                        "macierzy! (rows_substract)");
  assert(r2 < mat->r && "wiersz 2 (r2) jest jest wiekszy niz ilosc wierszow w "
                        "macierzy! (rows_substract)");

  for (int i = 0; i < mat->c; i++) {
    mat->data[r1][i] -= mat->data[r2][i];
  }
}

Matrix *create_extended_matrix(Matrix *A, Matrix *b) {
  mat_assert_square(A, __func__);
  mat_assert(b, __func__);
  assert(A->r == b->r && b->c == 1 &&
         "nie mozna rozszerzyc macierzy A o macierz b! Zle wymiary macierzy! "
         "(create_extended_matrix)");

  Matrix *extendedMatrix = createMatrix(A->r, A->c + b->c);
  mat_assert(extendedMatrix, __func__);

  for (int i = 0; i < A->r; i++) {
    for (int j = 0; j < A->c; j++) {
      extendedMatrix->data[i][j] = A->data[i][j];
    }
  }

  for (int i = 0; i < A->r; i++) {
    extendedMatrix->data[i][A->c] = b->data[i][0];
  }

  return extendedMatrix;
}

int find_max(Matrix *Ab, int currentColumn) {
  mat_assert(Ab, __func__);
  assert(currentColumn < Ab->c);

  int MaxRowIndex = 0;

  for (int i = 0; i < Ab->r; i++) {
    if (Ab->data[i][currentColumn] >= Ab->data[MaxRowIndex][currentColumn]) {
      MaxRowIndex = i;
    }
  }

  return MaxRowIndex;
}

int eliminate(Matrix *A, Matrix *b) {
  mat_assert_square(A, __func__);
  mat_assert(b, __func__);

  Matrix *Ab = create_extended_matrix(A, b);
  mat_assert(Ab, __func__);

  if (DEBUG)
    printToScreen(Ab);

  for (int pivot = 0; pivot < Ab->r; pivot++) {
    if (DEBUG)
      printf("Actual pivot: %lf \n", Ab->data[pivot][pivot]);

    for (int i = 1; i < Ab->r; i++) {
      if (pivot + i > Ab->r - 1) {
        break;
      }
      if (DEBUG)
        printf("Pivot+i:%d \n", pivot + i);

      if (!Ab->data[pivot][pivot]) {
        fprintf(stderr, "[!] Dzielenie przez 0! \n");
        return 1;
      }

      double scalar = (Ab->data[pivot + i][pivot] / Ab->data[pivot][pivot]);
      if (DEBUG)
        printf("Scalar: %lf \n", scalar);

      for (int j = 0; j < Ab->c; j++) {
        if (DEBUG)
          printf("Element %lf\n", Ab->data[pivot + i][j]);

        Ab->data[pivot + i][j] -= (scalar * Ab->data[pivot][j]);
      }
    }
  }
  if (DEBUG)
    printToScreen(Ab);

  if (DEBUG)
    printf("A->r:%d A->c:%d \n", A->r, A->c);
  for (int i = 0; i < A->r; i++) {
    for (int j = 0; j < A->c; j++) {
      A->data[i][j] = Ab->data[i][j];
      if (DEBUG)
        printf("i:%d j:%d \n", i, j);
    }
  }

  for (int i = 0; i < b->r; i++) {
    b->data[i][0] = Ab->data[i][Ab->c - 1];
  }

  freeMatrix(Ab);
  return 0;
}
double mod(double x){
  if(x<0) return -x;
  return x;
}

void swap_matrix(Matrix *mat, int c1, int c2) {
  mat_assert(mat, __func__);
  assert(c1 < mat->c && c2 < mat->c);

  Matrix *temp = createMatrix(mat->r, 1);
  mat_assert(temp, __func__);

  for (int i = 0; i < mat->r; i++) {
    temp->data[i][0] = mat->data[i][c1];
  }

  for (int i = 0; i < mat->r; i++) {
    mat->data[i][c1] = mat->data[i][c2];
  }

  for (int i = 0; i < mat->r; i++) {
    mat->data[i][c2] = temp->data[i][0];
  }

  if (DEBUG)
    printf("Swapped column nr. %d with column nr. %d\n", c1, c2);

  freeMatrix(temp);
}

int eliminate_B(Matrix *mat, Matrix *b) {
  mat_assert(b, __func__);
  mat_assert_square(mat, __func__);

  Matrix *Ab = create_extended_matrix(mat, b);
  mat_assert(Ab, __func__);

  if (DEBUG)
    printToScreen(Ab);

  for (int currentRow = 0; currentRow < Ab->r; currentRow++) {
    if (DEBUG)
      printf("Current Row: %d", currentRow);

    int MaxIndex = 0;
    for (int currentColumn = 0; currentColumn < Ab->c - 1; currentColumn++) {
      if (DEBUG)
        printf("Ab->data[%d][%d] = %lf \n", currentRow, currentColumn,
               Ab->data[currentRow][currentColumn]);
      if (mod(Ab->data[currentRow][currentColumn]) > mod(Ab->data[currentRow][MaxIndex]))
        MaxIndex = currentColumn;
    }
    if (DEBUG)
      printf("Ab->data[%d][MAX]: %lf, MAX=%d \n", currentRow,
             Ab->data[currentRow][MaxIndex], MaxIndex);

    if (MaxIndex != currentRow) {
      swap_matrix(Ab, currentRow, MaxIndex);
    }
    int pivot = currentRow;
    if (DEBUG)
      printf("Pivot: %d\n", pivot);

    if (DEBUG)
      printToScreen(Ab);

    for (int i = 1; i < Ab->r; i++) {
      if (Ab->data[pivot][pivot] == 0) {
        fprintf(stderr, "[!] Dzielenie przez 0! \n");
        return 1;
      }
      if (pivot + i >= Ab->r)
        break;
      double scalar = (Ab->data[pivot + i][pivot] / Ab->data[pivot][pivot]);

      if (DEBUG)
        printf("Scalar: %lf\n", scalar);

      for (int j = 0; j < Ab->c; j++) {
        if (DEBUG)
          printf("Element %lf\n", Ab->data[pivot + i][j]);
        Ab->data[pivot + i][j] -= (scalar * Ab->data[pivot][j]);
      }
    }
  }
  if (DEBUG)
    printf("A->r:%d A->c:%d \n", mat->r, mat->c);
  for (int i = 0; i < mat->r; i++) {
    for (int j = 0; j < mat->c; j++) {
      mat->data[i][j] = Ab->data[i][j];
      if (DEBUG)
        printf("i:%d j:%d \n", i, j);
    }
  }

  for (int i = 0; i < b->r; i++) {
    b->data[i][0] = Ab->data[i][Ab->c - 1];
  }
  if (DEBUG)
    printToScreen(Ab);

  freeMatrix(Ab);
  return 0;
}
