#ifndef _GAUSS_B_H
#define _GAUSS_B_H

#include "mat_io.h"

/**
 * Zwraca 0 - elimnacja zakonczona sukcesem
 * Zwraca 1 - macierz osobliwa - dzielenie przez 0
 */
int eliminate_B(Matrix *mat, Matrix *b);

#endif
