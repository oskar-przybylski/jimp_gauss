#include <stdio.h>
#include "backsubst.h"
/**
 * Zwraca 0 - wsteczne podstawienie zakonczone sukcesem
 * Zwraca 1 - błąd dzielenia przez 0 (element na diagonali = 0)
 * Zwraca 2 - błąd nieprawidłowych rozmiarów macierzy
 * x - macierz niewiadomych na której pracujemy
 * a - macierz wspolczynnikow
 * b - macierz wyrazow wolnych
 */
int  backsubst(Matrix *x, Matrix *a, Matrix *b) {
	for (int i = 0; i < x->r; i++) {
		x->data[i][0] = b->data[i][0];
	}


	double ratio;
	// int column = a->c - 1;
	// int column = 0;
	for(int column = a->c - 1; column >= 0; column--) {
		for(int i = a->r - (a->r - column) - 1; i >= 0; i--) {
			// a->data[i][column] = 0;

			//ratio
			a->data[i][column] *= a->data[column][column];
			ratio = a->data[i][column] / a->data[column][column];
			// a->data[i][column] -= ratio * a->data[column][column];
			printf("ratio: %lf\n", ratio);

			// for(int currentRow = i; currentRow >= 0; currentRow--) {
			// 	a->data[column][row]
			// }

			for(int currentCol = column - 1; currentCol >= 0; currentCol--) {
				a->data[i][currentCol] *= a->data[column][column];

			}

			a->data[i][column] -= ratio * a->data[column][column];


			printToScreen(a);
		}
	}

	printf("Zmodyfikowana macierz A:\n");
	printToScreen(a);


	return 0;
}


