#include "config.h"
#include "jacobi.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include <time.h>



int RunMethod(Jacobi *jacobi, j_type J_ERROR, int J_ITE_MAX);

void PrintUsage();

int main(int argc, char *argv[])
{
	Jacobi *jacobi = NULL;
	int i;
	int J_ORDER, J_ROW_TEST, J_ITE_MAX;
	j_type J_ERROR;
	int threads = 1;
	int error;
	FILE *fp = NULL;
	FILE *fpOut = NULL;
	int timesToExecute = 1;
	//! Execution times
	clock_t start;
	clock_t *end = NULL;
	double avarageTime = 0;
	double stdDeviation = 0;

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
	if( (end = (clock_t *) malloc(sizeof(clock_t) *timesToExecute)) == NULL )
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

	// Starting timer
	start = clock();

	// Running method
	for(i=0; i<timesToExecute; i++)
	{
		if( error = RunMethod(jacobi, J_ERROR, J_ITE_MAX) )
		{
			printf("Error #%dx!\n", error);
			goto end;
		}
		end[i] = clock();
	}
	
		
	// Calculate average and std. Deviation
	avarageTime = (double)(end[timesToExecute -1] - start) / CLOCKS_PER_SEC;
	avarageTime /= timesToExecute;

	stdDeviation = (double) (end[0] - start)/CLOCKS_PER_SEC -avarageTime;
	stdDeviation = stdDeviation *stdDeviation;
	for(i=1; i<timesToExecute; i++)
	{
		double tmp = (double) (end[i] - end[i-1])/CLOCKS_PER_SEC -avarageTime;
		stdDeviation += tmp *tmp;
	}
	stdDeviation = sqrt(stdDeviation);

	// Desired output
	printf("-------------------------------------------------------------------------------\n");
	printf("\t Average Execution Time:  %lf  +/- %lf \n", avarageTime, stdDeviation);
	printf("-------------------------------------------------------------------------------\n");
	printf("\t Iterations: %d \n", jacobi->iterations);
	printf("\t RowTest: %d => %lf =? %lf \n", J_ROW_TEST, jacobi->x1[J_ROW_TEST], jacobi->b[J_ROW_TEST]); //! TODO [fix this]
	printf("-------------------------------------------------------------------------------\n");

	// If user wants file with more output
	if(fpOut)
	{
		fprintf(fpOut, "Execution %2.d - Execution Time %lf \n", i, (double)(end[0] - start) / CLOCKS_PER_SEC);
		for(i=1; i<timesToExecute; i++)
			fprintf(fpOut, "Execution %2.d - Execution Time %lf \n", i, (double)(end[i] - end[i-1]) / CLOCKS_PER_SEC);
		fprintf(fpOut, "\t Average Execution Time:  %lf  +/- %lf \n", avarageTime, stdDeviation);
		fprintf(fpOut, "\t Iterations: %d \n", jacobi->iterations);
		fprintf(fpOut, "\t Solution: \n");
		for(i=0; i<jacobi->size; i++)
			fprintf(fpOut, "%lf\n", jacobi->b[i]);
	}

	// End
	end:
	free(end);
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
	printf("   ./prog <input file> [threads to use] [times to execute] [output file] \n");
	printf("      See README for input format \n");
	printf("      Default threads is 1 \n");
	printf("      Default time will execute is 1 \n");
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