#include "jacobi.h"
#include <stdlib.h>



int 	Jacobi_SinglePreprocess(j_type **A, j_type *b, int i, int size)
{
	j_type diagonal = A[i][i];
	int j;

	// Setting diagonal value to 0
	A[i][i] = 0;

	// Check if needed
	if(diagonal == 0)
		return 1;
	if(diagonal == 1)
		return 0;

	// Setting diagonal to 1
	for(j=0; j<size; j++)
		A[i][j] /= diagonal;
	b[i] /= diagonal;

	return 0;
}


j_type	Jacobi_SingleIteraction(j_type **A, j_type *b, j_type *x, int i, int size)
{
	j_type res = 0;
	int j;

	for(j=0; j<size; j++)
		res += A[i][j] *x[j];

	res = b[i] - res;

	return res;
}


j_type 	Jacobi_CheckPrecision(j_type *x1, j_type *x2, int size)
{
	int i;
	j_type precision;

	// Init precision
	precision = x2[0] - x1[0];
	precision = precision < 0 ? -precision : precision;

	// Calc largest deviation
	for(i=1; i<size; i++)
	{
		j_type tmp;

		tmp = x2[i] - x1[i];
		tmp = tmp < 0 ? -tmp : tmp;

		precision = precision < tmp ? tmp : precision;
	}

	return precision;
}