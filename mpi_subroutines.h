#include "parameters.h"
//#include <mpi.h>

/*-----require all subroutine header names to be included here------*/
void Analytical_Solution(double *x, double *y, double *Ttheor, double endTime); 
void Initialize_Solution(double *theta, double *Tnew, const int nRow, const int myRank, const int nProcs);
void update_Solution(double *theta, double *Tnew, const int nRow);

/*--------------MPI-routines-------------------------------------------------------------------------------*/
void mesh_2DGrid(double *x, double *y, int startRow, int nRow);
void countElements(const int N, const int nProcs, const int myRank, int *nRow, const int nGhostLayers);
void initData(double *theta, double *Ttheor, double *Tnew, const int myRank, const int start, const int end);
void exchange_SendRecv(double *theta,const int start, const int end, int src, int dest, const int myRank, const int nProcs);

void writeGlobalOutput(double *x, double *y, double *theta);
void writeProcData(double *x, double *y, double *theta, const int nRow, const int myRank, const int nProcs);
void writeProcMesh(double *x, double *y, const int nRow, const int myRank, const int nProcs);
void writeOutput(double *x, double *y, double *Ttheor, double *theta, const int myRank, const int start, const int end);
