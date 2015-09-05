#include "config.h"
#include "jacobi.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include <time.h>




void TimeToString(char *string, unsigned long long millisecs)
{
	unsigned int min = (unsigned int) (millisecs / (1000 *60) );
	millisecs -= (1000 *60) *min;
	unsigned int sec = (unsigned int) (millisecs / 1000);
	millisecs -= 1000 *sec;
	unsigned int msec = (unsigned int) millisecs;

	sprintf(string, "%.2u:%.2u,%.3u", min, sec, msec);
}



int RunMethod(Jacobi *jacobi, j_type J_ERROR, int J_ITE_MAX);

void PrintUsage();

int main(int argc, char *argv[])
{
	Jacobi *jacobi = NULL;
	int i;
	int J_ORDER, J_ROW_TEST, J_ITE_MAX;
	j_type J_ERROR;
	j_type rowTestVal = 0;
	int threads = 1;
	int error;
	FILE *fp = NULL;
	FILE *fpOut = NULL;
	int timesToExecute = 1;
	char tmp1[100], tmp2[100];
	//! Execution times
	struct timespec start;
	struct timespec end;
	unsigned long long int *execTime = NULL;
	unsigned long long int avarageTime;
	unsigned long long int stdDeviation;

	// Handling arguments
	if( argc < 2 )
	{
		PrintUsage();
		return 0;
	}
	if( argc >= 3 )
		threads = atoi(argv[2]);
	if( argc >= 4 )
		timesToExecute = atoi(argv[3]);
	if( argc >= 5 )
		if( (fpOut = fopen(argv[4], "w")) == NULL )
		{
			printf("Error opening output file!\n");
			return 0;
		};

	// Allocate all those execution times
	if( (execTime = (unsigned long long int *) malloc(sizeof(unsigned long long int) *timesToExecute)) == NULL )
	{
		printf("Error allocating memory!\n");
		goto end;
	}

	// Allocate internal structure
	if( (jacobi = (Jacobi *) malloc(sizeof(Jacobi))) == NULL )
	{
		printf("Error allocating memory!\n");
		goto end;
	}

	// Open file
	if( (fp = fopen(argv[1], "r")) == NULL )
	{
		printf("Error opening file!\n");
		goto end;
	}

	// Read initial info
	fscanf(fp, "%d", &J_ORDER);
	fscanf(fp, "%d", &J_ROW_TEST);
	fscanf(fp, "%lf", &J_ERROR);
	fscanf(fp, "%d", &J_ITE_MAX);

	// Init
	if( Jacobi_Init(jacobi, J_ORDER, threads) )
	{
		printf("Allocation error!\n");
		goto end;
	}

	// Read Matrices
	Jacobi_ReadMatrices(jacobi, fp);

	// Debug
	#ifdef DEBUG
	Jacobi_DebugMatrix(jacobi);
	Jacobi_DebugThreads(jacobi);
	#endif

	// Preprocess
	if( error = Jacobi_Preprocess(jacobi) )
	{
		printf("Error #%dp!\n", error);
		goto end;
	}


	// Running method
	avarageTime = 0;
	for(i=0; i<timesToExecute; i++)
	{
		// Starting timer
		clock_gettime(CLOCK_MONOTONIC, &start);

		// Run
		if( error = RunMethod(jacobi, J_ERROR, J_ITE_MAX) )
		{
			printf("Error #%dx!\n", error);
			goto end;
		}

		// Stop timer and calculate time
		clock_gettime(CLOCK_MONOTONIC, &end);
		execTime[i] = (unsigned long long int) 
		(
	 		(end.tv_sec - start.tv_sec) * 1000 + // Seconds to milliseconds 
			(end.tv_nsec - start.tv_sec) / 1000000 // Nanoseconds to milliseconds
		);
		avarageTime += execTime[i];
	}

	// Calculate average time
	avarageTime /= timesToExecute;

	// Calculate standard deviations
	stdDeviation = 0;
	for(i=0; i<timesToExecute; i++)
		stdDeviation += (execTime[i] - avarageTime)*(execTime[i] - avarageTime);
	stdDeviation = sqrt(stdDeviation);

	// Desired output
	#ifdef DIRTY_OUTPUT
	for(i=0; i<jacobi->size; i++)
		rowTestVal += jacobi->A[J_ROW_TEST][i] *jacobi->x1[i];
	rowTestVal += jacobi->x1[J_ROW_TEST]; // Cuz the diagonal was set to 0

	printf("-------------------------------------------------------------------------------\n");
	TimeToString(tmp1, avarageTime);
	TimeToString(tmp2, stdDeviation);
	printf("\t Average Execution Time: %s +/- %s \n", tmp1, tmp2);
	printf("-------------------------------------------------------------------------------\n");
	printf("\t Iterations: %d \n", jacobi->iterations);
	printf("\t RowTest: %d => %lf =? %lf \n", J_ROW_TEST, rowTestVal, jacobi->b[J_ROW_TEST]); //! TODO [fix this]
	printf("-------------------------------------------------------------------------------\n");
	#else
	TimeToString(tmp1, avarageTime);
	TimeToString(tmp2, stdDeviation);
	printf("\t Average Execution Time: %s +/- %s \n", tmp1, tmp2);
	printf("-------------------------------------------------------------------------------\n");
	#endif

	// If user wants file with more output
	if(fpOut)
	{
		for(i=1; i<timesToExecute; i++)
		{
			TimeToString(tmp1, execTime[i]);
			fprintf(fpOut, "Execution %2.d - Execution Time %s \n", i, tmp1);
		}
		TimeToString(tmp1, avarageTime);
		TimeToString(tmp2, stdDeviation);
		fprintf(fpOut, "\t Average Execution Time: %s +/- %s \n", tmp1, tmp2);
		fprintf(fpOut, "\t Iterations: %d \n", jacobi->iterations);
		fprintf(fpOut, "\t Solution: \n");
		for(i=0; i<jacobi->size; i++)
			fprintf(fpOut, "%lf\n", jacobi->b[i]);
	}

	// End
end:
	free(execTime);
	if(fpOut != NULL)
		fclose(fpOut);
	if(fp != NULL)
		fclose(fp);
	Jacobi_Destroy(jacobi);
	free(jacobi);
	return 0;
}


void PrintUsage()
{
	printf("Usage: \n");
	printf("   ./prog input-file-name [ number-of-threads [ number-of-times-to-execute [ output-file-name ] ] ] \n");
	printf("      See README for input format \n");
	printf("      Default threads is 1 \n");
	printf("      Default number of times will execute is 1 \n");
	printf("      Output will not be saved by default and contains more information \n");
}


int RunMethod(Jacobi *jacobi, j_type J_ERROR, int J_ITE_MAX)
{
	int i;

	// Init x1
	Jacobi_InitUnknowns(jacobi);
	jacobi->iterations = 0;

	// Debug
	#ifdef DEBUG
	Jacobi_Debug(jacobi);
	#endif

	// Jacobi Iterations
	return Jacobi_Run(jacobi, J_ERROR, J_ITE_MAX);
}