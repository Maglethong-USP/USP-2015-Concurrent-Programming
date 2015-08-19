#include "config.h"
#include "jacobi.h"

#include <stdio.h>
#include <stdlib.h>


#define size 3


int main(int argc, char *argv[])
{
	int i, j;
	j_type *A[size];
	j_type b[size];
	j_type x1[size];
	j_type x2[size];

	// alloc
	for(i=0; i<size; i++)
		A[i] = (j_type *) malloc(sizeof(j_type) * size);

	// Init
	A[0][0] = 10;	A[0][1] = 2;	A[0][2] = 1;	b[0] = 7;
	A[1][0] = 1;	A[1][1] = 5;	A[1][2] = 1;	b[1] = -8;
	A[2][0] = 2;	A[2][1] = 3;	A[2][2] = 10;	b[2] = 6;

	// Print matrix
	printf("Matrix: \n");
	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
			printf(" %8.4f", (float) (A[i][j]) );	
		printf("      %8.4f \n", (float) b[i]);
	}
	printf("\n");

	// Preprocess
	for(i=0; i<size; i++)
		if (Jacobi_SinglePreprocess(A, b, i, size))
		{
			printf("Error #001!\n");
			return 0;
		}

	// Print matrix
	printf("Matrix after processing: \n");
	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
			printf(" %8.4f", (float) (A[i][j]) );	
		printf("      %8.4f \n", (float) b[i]);
	}
	printf("\n");

	// Init x
	for(i=0; i<size; i++)
	{
		x1[i] = b[i];
		x2[i] = b[i];
	}

	// Print x
	printf("Iterations: \n[0] - ");
	for(i=0; i<size; i++)
		printf("%8.4f ", (float) x1[i] );
	printf("\n");



	for(j=0; j< MAX_ITER; j++)
	{
		j_type precision;

		// Iteration
		for(i=0; i<size; i++)
			x2[i] = Jacobi_SingleIteraction(A, b, x1, i, size);

		precision = Jacobi_CheckPrecision(x1, x2, size);

		// Swap
		for(i=0; i<size; i++)
			x1[i] = x2[i];

		// Print x
		printf("[%d] - ", j);
		for(i=0; i<size; i++)
			printf("%8.4f ", (float) x1[i] );
		printf(" - [%8.4f]", precision);
		printf("\n");

		// precision check
		if(precision < PRECISION)
			break;
	}


	for(i=0; i<size; i++)
		free(A[i]);
	return 0;
}