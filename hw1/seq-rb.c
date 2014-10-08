#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define MAX(a,b) ((a>b)? (a): (b))

double **grid;
int N, gridSize, MAXITERS;
double maxdiff;

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

void redblack() 
{
    int iters, i, j, jStart;
    double mydiff;

    for (iters = 1; iters <= MAXITERS+1; iters++) 
    {
        for (i = 1; i <= N; i++) 
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
	            maxdiff = MAX(maxdiff, fabs(grid[i][j] - mydiff));
	    }
    	}

	for (i = 1; i <= N; i++) 
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
	            maxdiff = MAX(maxdiff, fabs(grid[i][j] - mydiff));
	    }
    	}
    }
}


int main(int argc, char *argv[]) 
{
    int i, j;
    struct timeval tv;
    double startTime, endTime;

    if (argc != 3) 
    {
        exit(1);
    }

    gettimeofday(&tv, NULL);
    startTime = tv.tv_sec + tv.tv_usec/1000000.0;

    N = atoi(argv[1]);
    gridSize = N+2;

    MAXITERS = atoi(argv[2]);

    grid    = allocateGrid();

    /* Initialise grid including the boundaries */
    for (i = 0; i < gridSize; i++)
        for (j = 0; j < gridSize; j++) 
	{
	    if(i == 0 || i == gridSize-1 ||
	       j == 0 || j == gridSize-1)
		grid[i][j] = 1;
	    else grid[i][j] = 0;
    	}

    if (N <= 24)
        // print matrices if relatively small
        printGrid(grid);

    redblack();
    if (N <= 24)   // print matrix if relatively small
        printGrid(grid);
    
    gettimeofday(&tv, NULL);
    endTime = tv.tv_sec + tv.tv_usec/1000000.0;
    printf("#MPI Ranks : 0\t#Threads : 0\tExec. Time : %.3lf\tMaxdiff : %lf\n",
           (double)endTime - startTime, maxdiff);

}
