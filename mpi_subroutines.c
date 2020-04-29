#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> /* need it for MPI functions */
#include <sys/resource.h>
#include <stdbool.h>
#include <math.h>
#include "mpi_subroutines.h"
#include "timer.h"

/*-------subroutines for thermal 2d plate------------------*/
void mesh_2DGrid(double *x, double *y, int startRow, int nRow)
{
    int ic;
    for (int i = 0; i < nRow; i++){
        for (int j = 0; j < NCOL; j++) {

        ic = (i*NCOL) + j;
                                
        x[ic] = DX * (double) j;  // i-->rows BUT BE CAREFUL i increases in vertical (y) direction!!!
        y[ic] = DY * (double) (i+startRow);  // j-->columns BUT BE CAREFUL j increases in horizontal (x) direction!!!
        }
    }
}


void Analytical_Solution(double *x, double *y, double *Ttheor, double endTime) {

    int ic;
    double L = 1./LX;

	for (int i = 0; i < NROW; i++) {
	    for (int j = 0; j < NCOL; j++) {

               ic = (i*NCOL) + j;

	  	    double sum = 0.0;

                for (int n = 1; n < 101; n += 2) {
		    sum +=  (4.*T_N) / (n*PI) * sin(n*PI*x[ic]*L) * sinh(n*PI*y[ic]*L) / sinh(n*PI);
	        }
	        Ttheor[ic] = sum;
            }
	}
}



 
void Initialize_Solution(double *theta, double *Tnew, const int nRow, const int myRank, const int nProcs)
{
     int ic;
     // All other etries are correctly set at 0 by calloc
     // Set temperature at the internal nodes of the upper line of the rectangular domain
     if(myRank == nProcs -1)
     {
        int i = nRow-1;
        for (int j = 1; j < NCOL-1; j++){
               ic = (i*NCOL) + j;
               theta[ic]  = 100.0;
               Tnew  [ic] = 100.0;
        }

     }
}

void exchange_SendRecv(double *theta,const int start, const int end, int src, int dest, const int myRank, const int nProcs)
{
    int tag0 = 0;
    int tag1 = 1;
    int e = (end - 1) * NCOL;
    int s = 0;

    MPI_Sendrecv(&theta[e], NCOL, MPI_DOUBLE, dest, tag0, &theta[s], NCOL, MPI_DOUBLE, src, tag0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    s = (start + 1) * NCOL;
    e = end * NCOL;

    // NOTE WE ARE NOW SWAPPING THE SRC & DEST
    int tmp = dest;
    dest    = src;
    src     = tmp;

    MPI_Sendrecv(&theta[s], NCOL, MPI_DOUBLE, dest, tag1,&theta[e], NCOL, MPI_DOUBLE, src, tag1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}


void countElements(const int N, const int nProcs, const int myRank, int *nRow, const int nGhostLayers)
{
    int remainder = N % nProcs;
    if (remainder == 0) {
        *nRow = (N / nProcs) + nGhostLayers;
    }
    else {
        int pointsPerProcess = (N - remainder) / nProcs;
        if (myRank < remainder) {
            *nRow = pointsPerProcess + nGhostLayers + 1;
        }
        else
        {
           *nRow = pointsPerProcess + nGhostLayers;
        }
    }
}

void update_Solution(double *theta, double *Tnew, const int nRow)
{

    // Update internal cells
    int ic;
    double dd = 1.0 / (DXX);
    for (int i = 1; i < nRow-1; i++) {
        for (int j = 1; j < NCOL-1; j++){

           ic = (i*NCOL) + j;

           // compute future value for φ from eq(9) at time step t+1!
           Tnew[ic] = DT * dd * (theta[(i+1)*NCOL+j] + theta[(i-1)*NCOL+j] + theta[ic + 1] +\
                                           theta[ic - 1] -4. * theta[ic]) + theta[ic];
        }
    }
}

void writeOutput(double *x, double *y, double *Ttheor, double *theta, const int myRank, const int start, const int end)
{

    int ENOUGH = (int)(myRank/10+1+15);
    char fileName[ENOUGH];
    FILE *output;

    sprintf(fileName,"processor%d.mesh",myRank);
    output = fopen(fileName, "w");    // Φ(Χ,Υ) on 2D numerical mesh (41 * 21 = 861 nodes)

    

    int ic;
  //  int nrow = (end - start) + 1;
    for (int i = 0; i < NROW; i++){
        for (int j = 0; j < NCOL; j++) {

             ic = (i*NCOL) + j;
                     
           fprintf(output, "X = %10f     Y = %10f     T(theor) = %10f   theta(CPU) = %10f\n",\
                            x[ic],      y[ic],      Ttheor[ic],      theta[ic]  );
        }
    }

    fclose(output);
}


void writeProcMesh(double *x, double *y, const int nRow, const int myRank, const int nProcs)
{
    int ENOUGH = (int)(myRank/10+1+15);
    char fileName[ENOUGH];
    sprintf(fileName,"processor%d.mesh",myRank);

    printf("process %d: Writing mesh... in %s\n",myRank,fileName);
    FILE *output;
    output = fopen(fileName, "w");    // Φ(Χ,Υ) on 2D numerical mesh (41 * 21 = 861 nodes)

    int start=1;
    int end =nRow-1;
    if (myRank == 0 )
        start=0;
    if (myRank ==  nProcs -1 )
        end=nRow;

    int ic;
    for (int i = start; i < end; i++){
        for (int j = 0; j < NCOL; j++) {

             ic = (i*NCOL) + j;
                     
           fprintf(output, "X = %10f     Y = %10f\n",\
                            x[ic],      y[ic]);
        }
    }

    fclose(output);

}


void writeProcData(double *x, double *y, double *theta, const int nRow, const int myRank, const int nProcs)
{
    int ENOUGH = (int)(myRank/10+1+14);
    char fileName[ENOUGH];
    sprintf(fileName,"processor%d.dat",myRank);

    printf("process %d: Writing processor ... in %s\n",myRank,fileName);
    FILE *output;
    output = fopen(fileName, "w");    // Φ(Χ,Υ) on 2D numerical mesh (41 * 21 = 861 nodes)

    int start=1;
    int end =nRow-1;
    if (myRank == 0 )
        start=0;
    if (myRank ==  nProcs -1 )
        end=nRow;

    int ic;
    for (int i = start; i < end; i++){
        for (int j = 0; j < NCOL; j++) {

             ic = (i*NCOL) + j;
                     
           fprintf(output, "X = %10f     Y = %10f     theta(CPU) = %10f\n",\
                            x[ic],      y[ic],        theta[ic]);
        }
    }

    fclose(output);
}

void writeGlobalOutput(double *x, double *y, double *theta)
{
    FILE *output;
    output = fopen("Plate2DNumerical_mpi.dat", "w");    // Φ(Χ,Υ) on 2D numerical mesh (41 * 21 = 861 nodes)

    for (int i = 0; i < NROW*NCOL; i++){
     
       fprintf(output, "X = %10f     Y = %10f     theta(CPU) = %10f\n",\
                        x[i],        y[i],      theta[i]);
    }

    fclose(output);
}

