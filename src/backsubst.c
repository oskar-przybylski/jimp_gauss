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
	if(x == NULL || a == NULL || b == NULL) {
		fprintf(stderr, "[BŁĄD] Któraś z macierzy nie została zainicjalizowana!\n");
		return -1;
	}

	if(a->r != b->r) {
		fprintf(stderr, "[BŁĄD] Ilość rzędów macierzy współczynników nie jest równa ilości rzędów macierzy wyrazów wolnych!\n");
		return 2;
	}

	if(a->r < 1 || a->c < 1 || b->r < 1) {
		fprintf(stderr, "[BŁĄD] Nieprawidłowe wymiary macierzy!\n");
	}
	

	//Kopiowanie macierzy wyrazów wolnych b do x
	for (int i = 0; i < x->r; i++) {
		x->data[i][0] = b->data[i][0];
	}


	double ratio;
	for(int column = a->c - 1; column >= 0; column--) {
		for(int i = a->r - (a->r - column) - 1; i >= 0; i--) {
			//Zapobieganie dzieleniu przez 0
			if(a->data[column][column] == 0) {
				fprintf(stderr, "[BŁĄD] Dzielenie przez 0 - element na miejscu (%d, %d)!\n", column + 1, column + 1);
				return 1;
			}


			//Wymnożenie całego rzędu macierzy przez pivota - słabe podejście jeśli będą duże liczby, spróbuję potem poprawić
			for(int currentCol = column; currentCol >= 0; currentCol--) {
				a->data[i][currentCol] *= a->data[column][column];
			}

			//To samo tyle że dla macierzy x
			x->data[i][0] *= a->data[column][column];

			ratio = a->data[i][column] / a->data[column][column];

			//Odejmowanie rzędów w obu macierzach
			a->data[i][column] -= ratio * a->data[column][column];
			x->data[i][0] -= ratio * x->data[column][0];

			printToScreen(a);
			printToScreen(x);

		}
	}

	//Dzielenie macierzy wyników przez wartość danych współczynników z macierzy A
	for(int i = 0; i < x->r; i++) {
		if(a->data[i][i] == 0) {
			x->data[i][0] = 0;
			continue;
		}

		x->data[i][0] = x->data[i][0] / a->data[i][i];	
	}

	printf("Zmodyfikowana macierz A:\n");
	printToScreen(a);


	return 0;
}


