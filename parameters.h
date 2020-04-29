#ifndef PARAMS_H
#define PARAMS_H
/*----define constant parameters for the 2D thermal plate-----*/
#define LX 20.0 /* length in cm of the domain in x-direction  */
#define NCOL 101   /* includes boundary points on both end */
#define DX LX / (NCOL - 1.0)
#define LY 20.0 /* length in cm of the domain in y-direction  */
#define NROW 101   /* includes boundary points on both end */
#define DY LY / (NROW - 1.0)

#define DXX DX*DX   
#define DYY DY*DY

#define T_N 100.        // Top edge temperature maintained at 100 [C]
#define DT 0.225*DX*DY  //0.9*(1/4)*DX^2    

#define PI 4.*atan(1.) 
#endif
