#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TAG 13
#define MAX(a,b) ((a>b)? (a) : (b))

int main(int argc, char *argv[]) 
{
    double **grid, *tmp, maxdiff = 0.0, mydiff = 0.0, MAXDIFF = 0.0;
    double startTime, endTime;
    int numElements, offset, stripSize, gridSize, myrank; 
    int	HEIGHT, MAXITERS, numnodes, N, i, j, jStart, k;
    int firstRow, lastRow, iters;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    N = atoi(argv[1]);
    gridSize = N+2;
    MAXITERS = atoi(argv[2]);
    HEIGHT = N/numnodes;

    
    if (myrank == 0 && N<10) 
    {
    	tmp = (double *) malloc (sizeof(double ) * gridSize * gridSize);
        grid = (double **) malloc (sizeof(double *) * gridSize);
    	for (i = 0; i < gridSize; i++)
      	    grid[i] = &tmp[i * gridSize];
    }
    else 
    {
	tmp = (double *) malloc (sizeof(double ) * (HEIGHT + 2) * gridSize);
    	grid = (double **) malloc (sizeof(double *) * gridSize);
    	for (i = 0; i < gridSize; i++)
            grid[i] = &tmp[i * gridSize];
    } 
    
    firstRow = myrank * HEIGHT + 1;
    lastRow = firstRow + HEIGHT - 1;

    /* Initialise grid including the boundaries */
    for (i = 1; i <= HEIGHT; i++)
    {
        for (j = 1; j <= N; j++)
        {
	    grid[i][j] = 0;
            if(i == 1)
            {
                if(j == 1) grid[i-1][j-1] = 1;
                if(firstRow == 1) grid[i-1][j] = 1;
		    else grid[i-1][j] = 0;
            }
            if(i == HEIGHT)
            {
                if(j == N) grid[i+1][j+1] = 1;
                if(lastRow == N) grid[i+1][j] = 1;
		    else grid[i+1][j] = 0;
            }
            if(j == 1)
            {
                if(i == HEIGHT) grid[i+1][j-1] = 1;
                grid[i][j-1] = 1;
            }
            if(j == N)
            {
                if(i == 1) grid[i-1][j+1] = 1;
                grid[i][j+1] = 1;
            }
        }
    }

    /* Ensure that no node moves ahead until the entire grid is initialised */
    MPI_Barrier(MPI_COMM_WORLD);

    // start timer
    if (myrank == 0) 
    {
    	startTime = MPI_Wtime();
    }
    
    // do the work
    for (iters = 1; iters <= MAXITERS+1; iters++)
    {
        for (i = 1; i <= HEIGHT; i++)
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

	if(myrank > 0) 
	{
	    MPI_Send(grid[1], gridSize, MPI_DOUBLE, 
	             myrank-1, TAG, MPI_COMM_WORLD);
	    MPI_Recv(grid[0], gridSize, MPI_DOUBLE, 
		     myrank-1, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	} 
	if(myrank < numnodes-1)
	{
	    MPI_Recv(grid[HEIGHT+1], gridSize, MPI_DOUBLE,
		     myrank+1, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    MPI_Send(grid[HEIGHT], gridSize, MPI_DOUBLE, 
		     myrank+1, TAG, MPI_COMM_WORLD);
	}

	/* Sync the nodes to ensure symmetric values for the black computation */
	MPI_Barrier(MPI_COMM_WORLD);

	
        for (i = 1; i <= HEIGHT; i++)
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

	if(myrank > 0) 
	{
	    MPI_Send(grid[1], gridSize, MPI_DOUBLE, 
	             myrank-1, TAG, MPI_COMM_WORLD);
	    MPI_Recv(grid[0], gridSize, MPI_DOUBLE, 
		     myrank-1, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	} 
	if(myrank < numnodes-1)
	{
	    MPI_Recv(grid[HEIGHT+1], gridSize, MPI_DOUBLE,
		     myrank+1, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    MPI_Send(grid[HEIGHT], gridSize, MPI_DOUBLE, 
		     myrank+1, TAG, MPI_COMM_WORLD);
	}
	
	/* Sync for next iteration which begins with red computation */
	MPI_Barrier(MPI_COMM_WORLD);
    }
    /* Ensure all nodes reach this point before MPI_Reduce(maxdiff) is calculated in main() */
    MPI_Barrier(MPI_COMM_WORLD);

    // master receives from workers  -- note could be done via MPI_Gather
    if (N < 10 && myrank == 0) 
    {
    	offset = HEIGHT + 1;
        for (i=1; i<numnodes; i++) 
	{
	    if(i == numnodes-1) numElements = (HEIGHT+1) * gridSize;
	    	else numElements = HEIGHT * gridSize;
      	    
	    MPI_Recv(grid[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      	    offset += HEIGHT;
    	}
    }
    else if(N<10) 
    { 
    	// send my contribution to C
	if(myrank == numnodes - 1)
    	    MPI_Send(grid[1], (HEIGHT +1) * gridSize, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD);
    	MPI_Send(grid[1], HEIGHT * gridSize, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD);
    }

    MPI_Reduce(&maxdiff, &MAXDIFF, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // stop timer
    if (myrank == 0) 
    {
    	endTime = MPI_Wtime();
   	printf("#MPI Ranks : %d\t#Threads : 0\tExec. Time : %.3lf\tMaxdiff : %lf\n", numnodes,
           (double)endTime - startTime, MAXDIFF);
    }
    // print out matrix here, if I'm the master
    if (N < 10 && myrank == 0) 
    {
    	for (i=0; i < gridSize; i++) 
	{

		for (j=0; j < gridSize; j++) 
	    	{
        	    printf("%lf ", grid[i][j]);
      	    	}
		printf("\n");
    	}
    } 
    MPI_Finalize();
    return 0;
}

