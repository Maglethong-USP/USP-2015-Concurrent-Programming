/*
	Functions and structs starting with '_' are for internal use
*/
#ifndef __JACOBI__
#define __JACOBI__

#include <pthread.h>
#include <semaphore.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>		// for ceil()

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

	//! One mutex for each thread to sync to
	/*
		Each thread that calculates iterations will lock its own mutex
		after calculating one iteration going to sleep until the main
		thread unlocks it guaranteeing the thread will only calculate
		one iteration before the main thread orders another iteration.

		child 	locks it
		main 	unlocks it
	*/
	pthread_mutex_t *threadSync;
	//! One mutex so the threads don't write results at the same time
	/*
		Mainly to avoid cache bouncing when multiple threads try to
		write on the same array simultaneously. They all write in
		independent places, but they are all close to each other so
		the cache might go crazy and cause some delay.
	*/
	pthread_mutex_t resultSync;
	//! Semaphore for main thread to sync to children thread iterations
	/*
		The main thread will decrement the number of children threads
		from it while the child threads will increase it again
		guaranteeing all threads will calculate one iteration before
		the main thread continues.

		main decreases it [child count] times
		child increases once
	*/
	sem_t mainSync;
	
} Jacobi;

//! Struct with information for a thread
struct _Jacobi_ThreadInfo
{
	//! Reference to the Jacoby struct
	struct _Jacobi *jacobi;
	//! Start of the lines this thread will process
	int lineStart;
	//! End of the lines this thread will process
	int lineEnd;
	//! Number of this thread
	int num;
};



//! Verify if the method is guaranteed to converge
/*!
	/param jacobi			The initialized Jacobi struct with its matrices properly read
	/return					0 - Convergence guaranteed
							1 - Convergence not guaranteed
							2 - Diagonal as ZERO value
*/
int 	Jacobi_Verify(Jacobi *jacobi);


//! Initialize the Jacobi struct.
/*!
	This will not set the matrices, it will only prepare the struct.
	Use Jacobi_Destroy() to terminate it afterwards.
	/param jacobi			Reference to the Jacobi struct
	/param size				The size of the system
	/param threads			Number of threads that will run
	/return					0 		- No Error
							1 		- Memory allocation error
							> 10 	- Mutex initialization error
							> 100 	- Semaphore initialization error
*/
int 	Jacobi_Init(Jacobi *jacobi, int size, int threads);


//! Destroys the Jacoby struct and frees the used memory
/*!
	/param jacobi			The initialized Jacobi struct
*/
void 	Jacobi_Destroy(Jacobi *jacobi);


//! Prints debug output to screen containing thread information
/*!
	/param jacobi			The initialized Jacobi struct
*/
void 	Jacobi_DebugThreads(Jacobi *jacobi);


//! Prints debug output to screen containing information about both matrices and unknowns
/*!
	/param jacobi			The initialized Jacobi struct
*/
void 	Jacobi_Debug(Jacobi *jacobi);


//! Prints debug output to screen containing information about the constants and coefficients matrices
/*!
	/param jacobi			The initialized Jacobi struct
*/
void 	Jacobi_DebugMatrix(Jacobi *jacobi);


//! Prints debug output to screen containing information about the unknowns
/*!
	/param jacobi			The initialized Jacobi struct
*/
void 	Jacobi_DebugUnknowns(Jacobi *jacobi);


//! Initializes the unknowns with the current values of the constants matrix
/*!
	This is usually a good guess to start with the method
	/param jacobi			The initialized Jacobi struct
*/
void 	Jacobi_InitUnknowns(Jacobi *jacobi);


//! Processes the constants and the coefficients matrices fur further use
/*!
	Changes the matrix so that the diagonal values are 1.
	To further optimize the functions the diagonal values are set to 0 afterwards

		Ex: A = | 2 4 |	b = | 2 |	-> 	A = | 1 2 |	b = | 1 |	->	A = | 0 2 |	b =	| 1 |
				| 4 2 |		| 2 |			| 2 1 |		| 1 |			| 2 0 |		| 1 |

	/param jacobi			The initialized Jacobi struct
	/param return			0	if no error
							> 0	if error related to threads
*/
int 	Jacobi_Preprocess(Jacobi *jacobi);


//! Runs the Jacoby Method until the desired precision is reached
/*!
	/param jacobi			The initialized Jacobi struct
	/param desiredPrecision	The desired precision the method will try to reach
	/param maxIter			Maximum number of iterations that will run
	/param return			0	if no error
							> 0	if error related to threads
*/
int 	Jacobi_Run(Jacobi *jacobi, j_type desiredPrecision, int maxIter);


//! Job of one thread
/*!
	/param info				Thread information struct containing line range to process and reference to the Jacobi struct
*/
void 	_Jacobi_Thread(struct _Jacobi_ThreadInfo *info);


//! Processes some lines of the constants and the coefficients matrices fur further use
/*!
	/param info				Thread information struct containing line range to process and reference to the Jacobi struct
*/
void 	_Jacobi_ThreadPreprocess(struct _Jacobi_ThreadInfo *info);


//! [Internal] Calculates one iteration of one single unknown
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
j_type _Jacobi_SingleIteration(j_type **coefficients, j_type *constants, j_type *unknowns, int line, int size);


//! [Internal] Processes input matrices for other functions. Only one single line is preprocessed.
/*!
	Changes the matrix so that the diagonal values are 1.
	To further optimize the functions the diagonal values are set to 0 afterwards

		Ex: A = | 2 4 |	b = | 2 |	-> 	A = | 1 2 |	b = | 1 |	->	A = | 0 2 |	b =	| 1 |
				| 4 2 |		| 2 |			| 2 1 |		| 1 |			| 2 0 |		| 1 |

	/param coefficients		The matrix [A] containing the coefficients of the system.
	/param constants		The matrix [b] with the constant terms of the system.
	/param line				The line currently being preprocessed
	/param size				The size of the system (number of equations)
*/
void 	_Jacobi_SinglePreprocess(j_type **coefficients, j_type *constants, int line, int size);


//! [Internal] Checks the precision of the last iteration
/*!
	/param unknowns1 		The matrix [x] with the values of the unknown terms of the system of the previous.
	/param unknowns2 		The matrix [x] with the values of the unknown terms of the system of the last iteration.
	/param size				The size of the system (number of equations)
	/return					The precision of the last iteration's values
*/
j_type 	_Jacobi_CheckPrecision(j_type *unknowns1, j_type *unknowns2, int size);

#endif