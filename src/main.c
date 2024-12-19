#include "gauss_B.h"
#include "backsubst.h"
#include "mat_io.h"
#include <stdlib.h>
#include <stdio.h>


int main(int argc, char ** argv) {
	int res;
	Matrix * A = readFromFile(argv[1]);
	Matrix * b = readFromFile(argv[2]);
	Matrix * x;

	if (A == NULL) return -1;
	if (b == NULL) return -2;
  
  printf("Macierz A: \n");
	printToScreen(A);
  printf("Macierz b: \n");
	printToScreen(b);

	res = eliminate_B(A,b);
	x = createMatrix(b->r, 1);
	if (x != NULL) {
    
    if(res == 0){
		  res = backsubst(x,A,b);
    }else{
      fprintf(stderr,"Funkcja eliminate zwrocila wartosc inna niz 0! (Dzielenie przez 0 lub uklad sprzeczy) \n");
      return EXIT_FAILURE;
    } 
    
    if(res == 1){
       fprintf(stderr,"Funkcja backsubst zwrocila wartosc 1! (Dzielenie przez 0) \n");
       return EXIT_FAILURE;
    }else if(res == 2){
       fprintf(stderr,"Funkcja backsubst zwrocila wartosc 2! (Nieprawidlowe rozmiary macierzy) \n");
       return EXIT_FAILURE;
    }

    printf("Macierz x: \n");
		printToScreen(x);

	  freeMatrix(x);
	} else {
					fprintf(stderr,"Błąd! Nie mogłem utworzyć wektora wynikowego x.\n");
	}

	freeMatrix(A);
	freeMatrix(b);

	return 0;
}
