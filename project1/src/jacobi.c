#include "jacobi.h"


#include "config.h"


int 	Jacobi_Init(Jacobi *jacobi)
{
	int i;

	// Set all to NULL
	jacobi->iterations 	= 0;
	jacobi->threadSize 	= 0;
	jacobi->size 		= 0;
	jacobi->A 			= NULL;
	jacobi->b 			= NULL;
	jacobi->x1 			= NULL;
	jacobi->x2 			= NULL;
	jacobi->thread 		= NULL;
	jacobi->info 		= NULL;


	// Allocating

	// Threads
	jacobi->threadSize = THREADS;
	// Size
	jacobi->size = 3;
	// A
	if( (jacobi->A = (j_type **) malloc(sizeof(j_type *) *jacobi->size) ) == NULL )
		return 1;
	// Set A[i] NULL
	for(i=0; i<jacobi->size; i++)
		jacobi->A[i] = NULL;
	// A[i]
	for(i=0; i<jacobi->size; i++)
		if( (jacobi->A[i] = (j_type *) malloc(sizeof(j_type) *jacobi->size) ) == NULL )
			return 1;
	// b
	if( (jacobi->b = (j_type *) malloc(sizeof(j_type) *jacobi->size) ) == NULL )
		return 1;
	// x1
	if( (jacobi->x1 = (j_type *) malloc(sizeof(j_type) *jacobi->size) ) == NULL )
		return 1;
	// x2
	if( (jacobi->x2 = (j_type *) malloc(sizeof(j_type) *jacobi->size) ) == NULL )
		return 1;
	// thread
	if( (jacobi->thread = (pthread_t *) malloc(sizeof(pthread_t) *jacobi->threadSize) ) == NULL )
		return 1;
	// info
	if( (jacobi->info = (struct _Jacobi_ThreadInfo *) malloc(sizeof(struct _Jacobi_ThreadInfo) *jacobi->threadSize) ) == NULL )
		return 1;


	// Initializing

	// A
	jacobi->A[0][0] = 10;	jacobi->A[0][1] = 2;	jacobi->A[0][2] = 1;
	jacobi->A[1][0] = 1;	jacobi->A[1][1] = 5;	jacobi->A[1][2] = 1;
	jacobi->A[2][0] = 2;	jacobi->A[2][1] = 3;	jacobi->A[2][2] = 10;

	// b
	jacobi->b[0] = 7;
	jacobi->b[1] = -8;
	jacobi->b[2] = 6;

	// x1
	jacobi->x1[0] = 0;
	jacobi->x1[1] = 0;
	jacobi->x1[2] = 0;

	// x2
	jacobi->x2[0] = 0;
	jacobi->x2[1] = 0;
	jacobi->x2[2] = 0;

	// info
	for(i=0; i< jacobi->threadSize; i++)
	{
		jacobi->info[i].jacobi = jacobi;
		jacobi->info[i].line = 0;
	}

	// Success
	return 0; 
}

void 	Jacobi_Destroy(Jacobi *jacobi)
{
	int i;

	free(jacobi->info);
	free(jacobi->thread);
	free(jacobi->x2);
	free(jacobi->x1);
	free(jacobi->b);
	if(jacobi->A)
	{
		for(i=0; i<jacobi->size; i++)	
			free(jacobi->A[i]);
		free(jacobi->A);
	}
}

void 	Jacobi_Debug(Jacobi *jacobi)
{
	Jacobi_DebugMatrix(jacobi);
	Jacobi_DebugUnknowns(jacobi);
}

void 	Jacobi_DebugMatrix(Jacobi *jacobi)
{
	int i, j;

	// Print matrix
	printf("Matrix: \n");
	for(i=0; i<jacobi->size; i++)
	{
		for(j=0; j<jacobi->size; j++)
			printf(" %8.4f", (float) (jacobi->A[i][j]) );	
		printf("      %8.4f \n", (float) jacobi->b[i]);
	}
	printf("\n");
}

void 	Jacobi_DebugUnknowns(Jacobi *jacobi)
{
	int i;

	// Print x
	printf("Iteration %2.d - ", jacobi->iterations);
	for(i=0; i<jacobi->size; i++)
		printf("%8.4f ", (float) jacobi->x1[i] );
	printf("\n");
}

void 	Jacobi_InitUnknowns(Jacobi *jacobi)
{
	int i;
	for(i=0; i<jacobi->size; i++)
		jacobi->x1[i] = jacobi->b[i];
}

int 	Jacobi_Preprocess(Jacobi *jacobi)
{
	//! TODO [temporary only sequential]
	int i;
	struct _Jacobi_ThreadInfo info;

	// Init thread info
	info.jacobi = jacobi;

	// Start preprocessing each line
	for(i=0; i< jacobi->size; i++)
	{
		info.line = i;
		_Jacobi_SinglePreprocess2(&info);
	}

	return 0;
}

int 	Jacobi_Run(Jacobi *jacobi, j_type desiredPrecision)
{
	//! TODO [temporary only sequential]
	int i, j;
	j_type precision;
	struct _Jacobi_ThreadInfo info;
	j_type *tmp;

	// Init thread info
	info.jacobi = jacobi;

	// Run iterations
	for(j=0; j< MAX_ITER; j++)
	{
		// Iteration
		jacobi->iterations++;
		for(i=0; i< jacobi->size; i++)
		{
			info.line = i;
			_Jacobi_SingleIteraction2(&info);
		}
		// Precision calc
		precision = _Jacobi_CheckPrecision(jacobi->x1, jacobi->x2, jacobi->size);

		// Swap x1 with x2 [better than copying and erasing]
		tmp = jacobi->x1;
		jacobi->x1 = jacobi->x2;
		jacobi->x2 = tmp;

		// Debug
		Jacobi_DebugUnknowns(jacobi);

		// Precision reach check
		if(precision < desiredPrecision)
			break;
	}

	return 0;
}

int 	_Jacobi_SinglePreprocess2(struct _Jacobi_ThreadInfo *info)
{
	return _Jacobi_SinglePreprocess(info->jacobi->A, info->jacobi->b, info->line, info->jacobi->size);
}

int 	_Jacobi_SinglePreprocess(j_type **A, j_type *b, int i, int size)
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

void 	_Jacobi_SingleIteraction2(struct _Jacobi_ThreadInfo *info)
{
	info->jacobi->x2[info->line] = _Jacobi_SingleIteraction(info->jacobi->A, info->jacobi->b, info->jacobi->x1, info->line, info->jacobi->size);
}

j_type	_Jacobi_SingleIteraction(j_type **A, j_type *b, j_type *x, int i, int size)
{
	j_type res = 0;
	int j;

	for(j=0; j<size; j++)
		res += A[i][j] *x[j];

	res = b[i] - res;

	return res;
}


j_type 	_Jacobi_CheckPrecision(j_type *x1, j_type *x2, int size)
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