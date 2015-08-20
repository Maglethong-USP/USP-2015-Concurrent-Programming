/*
	Functions and structs starting with '_' are for internal use
*/
#ifndef __JACOBI__
#define __JACOBI__

//! Maximum number of iterations before error
#define MAX_ITER 20

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>

//! Defining type to be used on method
typedef double j_type;

//! Struct with system information
typedef struct _Jacobi
{
	//! The matrix [A] containing the coefficients of the system.
	j_type **A;
	//! The matrix [b] with the constant terms of the system.
	j_type *b;
	//! The matrix [x (k-1)] with the values of the unknown terms of the system of the previous iteration.
	j_type *x1;
	//! The matrix [x (k)] with the values of the unknown terms of the system of the last iteration.
	j_type *x2;
	//! The size of the system
	int size;
	//! The number of threads to be used
	int threadSize;
	//! Reference to each thread
	pthread_t *thread;
	//! Thread info [argument for thread]
	struct _Jacobi_ThreadInfo *info;
	//! Iteration counter
	int iterations;
} Jacobi;

//! Struct with information for thread
struct _Jacobi_ThreadInfo
{
	struct _Jacobi *jacobi;
	int line;
};



int 	Jacobi_Init(Jacobi *jacobi);
void 	Jacobi_Destroy(Jacobi *jacobi);
void 	Jacobi_Debug(Jacobi *jacobi);
void 	Jacobi_DebugMatrix(Jacobi *jacobi);
void 	Jacobi_DebugUnknowns(Jacobi *jacobi);
void 	Jacobi_InitUnknowns(Jacobi *jacobi);
int 	Jacobi_Preprocess(Jacobi *jacobi);
int 	Jacobi_Run(Jacobi *jacobi, j_type desiredPrecision);
int 	_Jacobi_SinglePreprocess2(struct _Jacobi_ThreadInfo *info);
void 	_Jacobi_SingleIteraction2(struct _Jacobi_ThreadInfo *info);

//! Calculates one iteration of one single unknown
/*!
	/param coefficients		The matrix [A] containing the coefficients of the system.
								This matrix must have been pre-processed using Jacobi_SinglePreprocess().
	/param constants		The matrix [b] with the constant terms of the system.
								This matrix must have been pre-processed using Jacobi_SinglePreprocess().
	/param unknowns			The matrix [x] with the values of the unknown terms of the system of the last iteration.
								For the first iteration this matrix must have been pre-initialized with guessed values.
	/param line				The line of the unknown currently being updated
	/param size				The size of the system (number of equations)
	/return					The updated value of the unknown of this line
*/
j_type _Jacobi_SingleIteraction(j_type **coefficients, j_type *constants, j_type *unknowns, int line, int size);


//! Reprocess input matrices for other functions. Only one single line is preprocessed.
/*!
	Changes the matrix so that the diagonal values are 1.
	To further optimize the functions the diagonal values are set to 0 afterwards

		Ex: A = | 2 4 |	b = | 2 |	-> 	A = | 1 2 |	b = | 1 |	->	A = | 0 2 |	b =	| 1 |
				| 4 2 |		| 2 |			| 2 1 |		| 1 |			| 2 0 |		| 1 |

	/param coefficients		The matrix [A] containing the coefficients of the system.
	/param constants		The matrix [b] with the constant terms of the system.
	/param line				The line currently being preprocessed
	/param size				The size of the system (number of equations)
	/return					Error value: 
							 * 0 -> no error
							 * 1 -> Diagonal value is 0 (zero div error)
*/
int 	_Jacobi_SinglePreprocess(j_type **coefficients, j_type *constants, int line, int size);


//! Checks the precision of the last iteration
/*!
	/param unknowns1 		The matrix [x] with the values of the unknown terms of the system of the previous.
	/param unknowns2 		The matrix [x] with the values of the unknown terms of the system of the last iteration.
	/param size				The size of the system (number of equations)
	/return					The precision of the last iteration's values
*/
j_type 	_Jacobi_CheckPrecision(j_type *unknowns1, j_type *unknowns2, int size);

#endif