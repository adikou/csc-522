#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>

#define MAX(a,b) ((a>b)? (a): (b))

double **grid, *maxdiff;
int    N, HEIGHT, gridSize, MAXITERS, numThreads;
int    *arrive = 0;

double **allocateGrid()
{
    int i;
    double *vals, **temp;

    // allocate values
    vals = (double *) malloc (gridSize * gridSize * sizeof(double));

    // allocate vector of pointers
    temp = (double **) malloc (gridSize * sizeof(double*));

    for(i=0; i < gridSize; i++)
        temp[i] = &(vals[i * gridSize]);

    return temp;
}

void printGrid(double **grid)
{
    int i,j;

    printf("\nThe %d * %d grid is\n", gridSize, gridSize);
    for(i=0; i < gridSize; i++)
    {
        for(j=0; j < gridSize; j++)
            printf("%lf ",  grid[i][j]);
        printf("\n");
    }
}

void barrier(int id)
{
    int j, lookAt;
    for(j = 1; j <= ceil(log(numThreads)); j++)
    {	
	while(arrive[id] != 0);
	arrive[id] = j;
	lookAt = (int)(id + pow(2, j-1)) % numThreads;
	while(arrive[lookAt] != j);
	arrive[lookAt] = 0;
    }
    
}

void redblack(int id)
{
    int iters, i, j, firstRow, lastRow, jStart;
    double mydiff;

    firstRow = id * HEIGHT + 1;
    lastRow = firstRow + HEIGHT - 1;

    /* Initialise grid including the boundaries */
    for (i = firstRow; i <= lastRow; i++)
    {
        for (j = 1; j <= N; j++)
        {
            if(i == 1) 
	    {
		if(j == 1) grid[i-1][j-1] = 1;
		grid[i-1][j] = 1; 
	    }
	    else if(i == N) 
	    {
		if(j == N) grid[i+1][j+1] = 1;
		grid[i+1][j] = 1; 
	    }
	    else grid[i][j] = 0;
            if(j == 1) 
	    {
		if(i == N) grid[i+1][j-1] = 1;
		grid[i][j-1] = 1; 
	    }
	    else if(j == N) 
	    {
		if(i == 1) grid[i-1][j+1] = 1;
		grid[i][j+1] = 1;
	    }
            else grid[i][j] = 0;
        }
    }
 
    /* Ensure that no thread moves ahead until the entire grid is initialised */
    barrier(id);

    for (iters = 1; iters <= MAXITERS+1; iters++)
    {
        for (i = firstRow; i <= lastRow; i++)
        {
            if(i%2 == 1)  jStart = 1;
                else jStart = 2;


            for (j = jStart; j <= N; j += 2)
            {
                if(iters == MAXITERS + 1)
                     mydiff = grid[i][j];

                grid[i][j] = (grid[i-1][j] + grid[i][j-1] +
                             grid[i+1][j] + grid[i][j+1]) * 0.25;

                if(iters == MAXITERS + 1)
                    maxdiff[id] = MAX(maxdiff[id], fabs(grid[i][j] - mydiff));
            }
        }

	/* Sync the threads to ensure symmetric values for the black computation */
	barrier(id);
        for (i = firstRow; i <= lastRow; i++)
        {
            if(i%2 == 1)  jStart = 2;
                else jStart = 1;

            for (j = jStart; j <= N; j += 2)
            {
                if(iters == MAXITERS + 1)
                   mydiff = grid[i][j];

                grid[i][j] = (grid[i-1][j] + grid[i][j-1] +
                             grid[i+1][j] + grid[i][j+1]) * 0.25;

                if(iters == MAXITERS + 1)
                    maxdiff[id] = MAX(maxdiff[id], fabs(grid[i][j] - mydiff));
            }
        }

	/* Sync for next iteration which begins with red computation */
	barrier(id);
    }

    /* Ensure all threads reach this point before max(maxdiff) is calculated in main() */
    barrier(id);
}

void *worker(void *arg)
{
    int id = *((int *) arg);
    redblack(id);
    return NULL;
}

int main(int argc, char *argv[])
{
    int i, j;
    int *p;
    pthread_t *threads;
    double MAXDIFF = 0;
    struct timeval tv;
    double startTime, endTime;

    if (argc != 4)
    {
	printf("Usage: %s <size> <MAXITERS> <numThreads>,  where size" 
		"is dimension of grid matrix, MAXITERS is max iterations" 
		"and n is number of threads\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    gridSize = N+2;
    MAXITERS = atoi(argv[2]);
    numThreads = atoi(argv[3]);

    HEIGHT = N/numThreads;
    grid    = allocateGrid();

    maxdiff = (double*) malloc(numThreads * sizeof(double));
    arrive  = (int*)malloc(numThreads * sizeof(int));

    /* Initialise arrive and maxdiff arrays */
    for(i = 0; i < numThreads; i++)
    {
	arrive[i] = 0;
	maxdiff[i] = 0.0;
    }

    // Allocate thread handles
    threads = (pthread_t *) malloc(numThreads * sizeof(pthread_t));

    gettimeofday(&tv, NULL);
    startTime = tv.tv_sec + tv.tv_usec/1000000.0;

    // Create threads 
    for (i = 0; i < numThreads; i++) 
    {
    	p = (int *) malloc(sizeof(int));  // yes, memory leak, don't worry for now
    	*p = i;
    	pthread_create(&threads[i], NULL, worker, (void *)(p));
    }

    for (i = 0; i < numThreads; i++) 
    {
    	pthread_join(threads[i], NULL);
    }
    for (i = 0; i < numThreads; i++)
	MAXDIFF = MAX(MAXDIFF, maxdiff[i]);
    gettimeofday(&tv, NULL);
    endTime = tv.tv_sec + tv.tv_usec/1000000.0;

    if(N <= 10)
    	printGrid(grid);

    printf("#MPI Ranks : 0\t#Threads : %d\tExec. Time : %.3lf\tMaxdiff : %lf\n", numThreads, 
	   (double)endTime - startTime, MAXDIFF);

}
