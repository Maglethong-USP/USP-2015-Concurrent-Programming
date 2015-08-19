#ifndef __JACOBI__
#define __JACOBI__


//! Defining type to be used on method
typedef double j_type;

//! Calculates one iteraction of one single unknown
/*!
	/param coefficients		The matrix [A] containing the coefficients of the system.
								This matrix must have been pre-processed using Jacobi_SinglePreprocess().
	/param constants		The matrix [b] with the constant terms of the system.
								This matrix must have been pre-processed using Jacobi_SinglePreprocess().
	/param unknowns			The matrix [x] with the values of the unknown terms of the system of the last iteraction.
								For the first iteraction this matrix must have been pre-initialized with guessed values.
	/param line				The line of the unknown currently being updated
	/param size				The size of the system (number of equations)
	/return					The updated value of the unknown of this line
*/
j_type Jacobi_SingleIteraction(j_type **coefficients, j_type *constants, j_type *unknowns, int line, int size);


//! Preprocesses input matrixes for other functios. Only one single line is preprocessed.
/*!
	Changes the matrix so that the diagonal values are 1.
	To further optimize the functions the diagonal values are set to 0 afterwars

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
int 	Jacobi_SinglePreprocess(j_type **coefficients, j_type *constants, int line, int size);


//! Checks the precision of the last iteraction
/*!
	/param unknowns1 		The matrix [x] with the values of the unknown terms of the system of the previous.
	/param unknowns2 		The matrix [x] with the values of the unknown terms of the system of the last iteraction.
	/param size				The size of the system (number of equations)
	/return					The precision of the last iteraction's values
*/
j_type 	Jacobi_CheckPrecision(j_type *unknowns1, j_type *unknowns2, int size);

#endif