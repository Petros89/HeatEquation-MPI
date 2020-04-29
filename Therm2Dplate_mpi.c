/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * PROJECT 3 *           DUE DATE:: DECEMBER 14                      *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *           PARALLEL VERSION OF 2D HEAT EQUATION (MPI)              *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  Numerical & Analytical Heat Diffusion in a 2D Square Plate of    *
 *  side size = 20cm and north edge temperature TN = 100 C.          *
 *  Thermal conductivity is assumed k=1 for this problem             * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  Author: Petros Apostolou                                         *
 *  Date:   12/13/2018                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *  compile: make                                                    *
 *  execute: ./Heat2D.exe < simulation end time (seconds)>           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> /* need it for MPI functions */
#include <stdbool.h>
#include <math.h>
#include <sys/resource.h>
#include "parameters.h"
#include "mpi_subroutines.h"
#include "timer.h"


/*---main program for the parallel implementation of the 2D Heat Plate Problem*/
int main(int argc, char *argv[])
{
    if (argc < 2) {
        perror("Command-line usage: executableName < End Time (sec)>");
        exit(1);
    }
    
    double endTime = (double) atof(argv[ 1 ]);
    double time = 0.;
    double *temp;
    

    double *x ;
    double *y ;
    double *theta  ;
    double *Tnew   ;
    double *Ttheor ;
      

/*-------starting the parallelization process-----------------------------*/
    int nProcs; /* number of processes */
    int myRank; /* process rank */
    int src;    /* handles for communication, source process id */
    int dest;   /* handles for communication, destination process id */
    int startRow = 0 ;  /* start Row index for each partial domain */
    int endRow;    /* end Row index for each partial domain */
    int nRow;    /* end index for each partial domain */

    MPI_Init(&argc, &argv);                 /* initialize MPI */
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs); /* get the number of processes */

    int *procRowCount = (int*) calloc(nProcs, sizeof(int)); //pointer to the rows of each processor
    int nDims = 1;          // dimension of Cartesian decomposition: here 1 slice (horizontal slicing)
    int dimension[nDims];   // dimension of slicing
    int isPeriodic[nDims];  // dimension of periodicity in the domain
    int reorder = 1;        // allow system to optimize(reorder) the mapping of processes to physical cores

    dimension[0]  = nProcs; 
    isPeriodic[0] = 0;  // periodicty of each dimension

    MPI_Comm comm1D; // define a communicator that would be assigned a new topology
    MPI_Cart_create(MPI_COMM_WORLD, nDims, dimension, isPeriodic, reorder, &comm1D);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); /* get the rank of a process after REORDERING! */
    MPI_Cart_shift(comm1D, 0, 1, &src, &dest); /* Let MPI find out the rank of processes for source and destination */


/* let's make sure to insert adequate number of ghost cells for data communication between two consecutive processors */
    int nGhostLayers;
    if (myRank == 0 || myRank == nProcs - 1) {
        nGhostLayers = 1;
    } else {
        nGhostLayers = 2;
    }


/* this function counts the number of elements are treated by each processor */
    countElements(NROW, nProcs, myRank, &nRow, nGhostLayers);
    printf("my rank: %d, nRow: %d \n", myRank,nRow);


/* the command below gathers all data from all processes and distrubutes the combined data to all ranks */
    MPI_Allgather(&nRow, 1, MPI_INT, procRowCount,1,MPI_INT,MPI_COMM_WORLD);
    if (myRank == 0)
    {
        for (int i =0; i < nProcs; ++i)
        {

            printf("processor: %d: %d\n",i,procRowCount[i]);
        }
    }

    if (myRank > 0)
    {
        for (int i =0; i<myRank; ++i)
        {
            startRow += procRowCount[i]-2;
        }
    }

    endRow = startRow + nRow-1;


/* let's make sure that each processor takes care of each own range of rows [startRow, endRow] */
    printf("my rank: %d, nRow: %d, startRow: %d, endRow: %d\n", myRank,nRow,startRow,endRow);


    //Allocate resources
    x = (double*) calloc (nRow*NCOL, sizeof (*x));
    y = (double*) calloc (nRow*NCOL, sizeof (*y));
    theta  = (double*) calloc (nRow*NCOL, sizeof (*theta));
    Tnew   = (double*) calloc (nRow*NCOL, sizeof (*Tnew));
    Ttheor = (double*) calloc (nRow*NCOL, sizeof (*Ttheor));


    // calculate the  coordinates of each computational point
    mesh_2DGrid(x, y, startRow, nRow);


    //print output to .test file to see if each processor's mesh is correct
    writeProcMesh(x,y,nRow, myRank, nProcs);


    // Initialize Numerical Solution 
    Initialize_Solution(theta, Tnew,nRow, myRank, nProcs);


    // number of iterations
    int iters=0;


    // make sure all processors are ready to start the processing
    MPI_Barrier(MPI_COMM_WORLD); // barrier to start counting time


    // Time to start timing
    double MPI_start = MPI_Wtime();


/* update in the numerical solution starts below */
    while (time < endTime) { 
        
        // Send and Recieve datas simultaneously
        exchange_SendRecv(theta, 0, nRow-1, src, dest, myRank, nProcs);

        // update Temperature values
        update_Solution(theta, Tnew, nRow);
        
        temp  = Tnew;
        Tnew  = theta;
        theta = temp;

        time += DT;
        ++iters;

    }

    
    MPI_Barrier(MPI_COMM_WORLD); // barrier to get the finish time
    double MPI_finish = MPI_Wtime();

    
    writeProcData(x,y,theta,nRow, myRank, nProcs);

    // allocate global variables
    double *globalT;
    double *globalX;
    double *globalY;


    if (myRank == 0)
    {
        globalT = (double*)calloc(NROW*NCOL, sizeof(*globalT));
        globalX = (double*)calloc(NROW*NCOL, sizeof(*globalX));
        globalY = (double*)calloc(NROW*NCOL, sizeof(*globalY));
    }

    int startLocal = 1;
    if (myRank == 0)
        startLocal = 0;

    int *recvCounts = (int*) calloc(nProcs, sizeof(int));
    int sentCount   = (nRow-nGhostLayers)*NCOL;


    MPI_Allgather(&sentCount, 1, MPI_INT, recvCounts,1,MPI_INT,MPI_COMM_WORLD);

    int *displ = (int*) calloc(nProcs, sizeof(int));

    for (int i = 1; i < nProcs; i++)
    {
        displ[i] = displ[i-1] + recvCounts[i-1];
    }

/*--let's gather all processes to start printing results at the same time--*/
    MPI_Gatherv(&theta[startLocal*NCOL], sentCount, MPI_DOUBLE, globalT,
                recvCounts, displ, MPI_DOUBLE ,0,MPI_COMM_WORLD);
    MPI_Gatherv(&x[startLocal*NCOL],sentCount,MPI_DOUBLE, globalX,
                recvCounts, displ, MPI_DOUBLE ,0,MPI_COMM_WORLD);
    MPI_Gatherv(&y[startLocal*NCOL],sentCount,MPI_DOUBLE, globalY,
                recvCounts, displ, MPI_DOUBLE ,0,MPI_COMM_WORLD);

    if (myRank==0)

/* let's use first processor to write the globsl results in output .dat file */      
        writeGlobalOutput(globalX, globalY, globalT);


    double MPI_elapsedTime = MPI_finish - MPI_start;
    double wallTime;  // returns time in seconds !Be careful to print in ms (*1000.)


/*--let's make sure to time the maximum process time--*/
    MPI_Reduce(&MPI_elapsedTime, &wallTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myRank == 0){
        printf("Wall-clock time = %4.4f ms \n", wallTime*1000.);
        //printf("Executed %d iterations in total and simulated time is %lf ms \n", iters,time);
    }


/* let's free out the memory */
    free(x);
    free(y);
    free(theta);
    free(Tnew);
    free(Ttheor);

/* MPI's functionality ends here */
    MPI_Finalize( );

 return EXIT_SUCCESS;
}
