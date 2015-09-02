#include "jacobi.h"


#include "config.h"


int 	Jacobi_ReadMatrices(Jacobi *jacobi, FILE *fp)
{
	int error = 0;
	int i, j;


	// Read Matrix A
	for(i=0; i<jacobi->size; i++)
		for(j=0; j<jacobi->size; j++)
			fscanf(fp, "%lf", &(jacobi->A[i][j]));

	// Read Matrix b
	for(i=0; i<jacobi->size; i++)
		fscanf(fp, "%lf", &(jacobi->b[i]));

	// End
	return error;
}

int 	Jacobi_Init(Jacobi *jacobi, int size, int threads)
{
	int i;
	int linesPerThread;
	int error;

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
	jacobi->threadSync 	= NULL;


	// Allocating

	// Threads
	jacobi->threadSize 	= threads;
	// Size
	jacobi->size 		= size;

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

	// info
	linesPerThread = ceil( (float)jacobi->size / (float)jacobi->threadSize);
	for(i=0; i< jacobi->threadSize; i++)
	{
		jacobi->info[i].jacobi = jacobi;
		jacobi->info[i].lineStart = i *linesPerThread;
		jacobi->info[i].lineEnd = (i +1) *linesPerThread -1;
		jacobi->info[i].num = i;
	}
	jacobi->info[jacobi->threadSize -1].lineEnd = jacobi->size -1;

	// Mutexes
	if( (jacobi->threadSync = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t) *jacobi->threadSize) ) == NULL )
		return 1;

	for(i=0; i<jacobi->threadSize; i++)
	{
		if( error = pthread_mutex_init(jacobi->threadSync +i, NULL) )
			return error +10;
		// Lock mutex [Want my thread to start asleep]
		pthread_mutex_lock(jacobi->threadSync +i);
	}

	if( error = pthread_mutex_init( &(jacobi->resultSync), NULL) )
		return error +10;

	// Semaphore
		if( error = sem_init( &(jacobi->mainSync), 0, 0 ) )
			return error +100;

	// Success
	return 0; 
}

int 	Jacobi_Verify(Jacobi *jacobi)
{
	int i, j;
	int diverge = 0;
	int sum;

	for(i=0; i<jacobi->size; i++)
	{
		// 0 value in diagonal
		if(jacobi->A[i][i] == 0)
			diverge = 2;
		// Diagonal value must be superior to the rest of the line added up
		sum = 0;
		for(j=0; j<jacobi->size; j++)
			sum += jacobi->A[i][j];
		sum -= jacobi->A[i][i] *2;
		if(sum > 0)
			diverge = 1;
	}

	// 0 -> Converges | 1 May Diverge | 2 ZERO div
	return diverge;
}

void 	Jacobi_Destroy(Jacobi *jacobi)
{
	int i;

	
	sem_destroy( &(jacobi->mainSync) );
	pthread_mutex_destroy( &(jacobi->resultSync) );
	if(jacobi->threadSync)
	{
		for(i=0; i<jacobi->threadSize; i++)	
			pthread_mutex_destroy(jacobi->threadSync +1);
		free(jacobi->threadSync);
	}
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
	jacobi->threadSync 	= NULL;
}

void 	Jacobi_DebugThreads(Jacobi *jacobi)
{
	int i;

	printf("Threads: (%d) [%d lines]\n", jacobi->threadSize, jacobi->size);
	for(i=0; i<jacobi->threadSize; i++)
	{
		printf(" %2.d - [%d, %d] \n", i, jacobi->info[i].lineStart, jacobi->info[i].lineEnd);
	}
	printf("\n");
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
	jacobi->iterations = 0;
}

int 	Jacobi_Preprocess(Jacobi *jacobi)
{
	int i;
	int error = 0;

	// Create threads
	for(i=0; i< jacobi->threadSize && !error; i++)
		error = pthread_create(jacobi->thread +i, NULL, (void * (*)(void *)) _Jacobi_ThreadPreprocess, (void *)(jacobi->info +i));

	// Join threads to synchronize
	for(i=0; i< jacobi->threadSize && !error; i++)
		error = pthread_join(jacobi->thread[i], NULL);

	return error;
}

int 	Jacobi_Run(Jacobi *jacobi, j_type desiredPrecision, int maxIter)
{
	int i, j;
	j_type precision;
	j_type *tmp;
	int error = 0;

	// Create threads
	for(i=0; i< jacobi->threadSize && !error; i++)
		// Create thread
		error = pthread_create(jacobi->thread +i, NULL, (void * (*)(void *)) _Jacobi_Thread, (void *)(jacobi->info +i));

	// Run iterations
	for(i=0; i< maxIter; i++)
	{
		jacobi->iterations++;

		// Unlock each thread's mutex allowing it a new iteration
		for(j=0; j<jacobi->threadSize; j++)
			pthread_mutex_unlock(jacobi->threadSync +j);

		// Lower semaphore [child thread count] times to wait till they all finish
		for(j=0; j<jacobi->threadSize; j++)
			sem_wait( &(jacobi->mainSync) );

		// Precision calc
		precision = _Jacobi_CheckPrecision(jacobi->x1, jacobi->x2, jacobi->size);

		// Swap x1 with x2 [better than copying and erasing]
		tmp = jacobi->x1;
		jacobi->x1 = jacobi->x2;
		jacobi->x2 = tmp;

		// Debug
		#ifdef DEBUG
		Jacobi_DebugUnknowns(jacobi);
		#endif

		// Precision reach check
		if(precision < desiredPrecision || error)
			break;
	}

	// Kill all child threads
	for(i=0; i<jacobi->threadSize; i++)
		pthread_cancel(jacobi->thread[i]);

	return error;
}

void 	_Jacobi_Thread(struct _Jacobi_ThreadInfo *info)
{
	int i, j;
	j_type *tmp;	// To temporarily store results avoiding cache bouncing

	// Allocating
	if( (tmp = (j_type *) malloc( sizeof(j_type) *(info->lineEnd - info->lineStart +1) )) == NULL )
		return;

	// Run indefinitely. Main will have to kill me
	while(1)
	{
		// Lock this thread's mutex allowing only one iteration
		pthread_mutex_lock(info->jacobi->threadSync +(info->num));

		// Calculate one iteration
		//_Jacobi_ThreadIteration(info);

		for(i = info->lineStart, j=0; i <= info->lineEnd; i++, j++)
			tmp[j] = _Jacobi_SingleIteration(info->jacobi->A, info->jacobi->b, info->jacobi->x1, i, info->jacobi->size);

		// Copy temporary array to results
		pthread_mutex_lock( &(info->jacobi->resultSync) );
		for(i = info->lineStart, j=0; i <= info->lineEnd; i++, j++)
			info->jacobi->x2[i] = tmp[j];
		pthread_mutex_unlock( &(info->jacobi->resultSync) );

		// Raise the semaphore
		sem_post( &(info->jacobi->mainSync) );
	}

	// Free
	free(tmp);
}

void 	_Jacobi_ThreadPreprocess(struct _Jacobi_ThreadInfo *info)
{
	int i;

	for(i = info->lineStart; i <= info->lineEnd; i++)
		_Jacobi_SinglePreprocess(info->jacobi->A, info->jacobi->b, i, info->jacobi->size);
}

void 	_Jacobi_SinglePreprocess(j_type **A, j_type *b, int i, int size)
{
	j_type diagonal = A[i][i];
	int j;

	// Setting diagonal value to 0
	A[i][i] = 0;

	// Setting diagonal to 1
	for(j=0; j<size; j++)
		A[i][j] /= diagonal;
	b[i] /= diagonal;
}

j_type	_Jacobi_SingleIteration(j_type **A, j_type *b, j_type *x, int i, int size)
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