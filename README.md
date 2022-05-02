# Inter-node MPI Implementation of 3D Heat Transfer

## Uaage
- Load modules: module load gcc/5.4.0 ==> loads the GNU compiler for the gcc/5.4.0 version. module load mpich/3.1 ==> loads the mpich/3.1 distribution
- Compile: make  ==> links all object files and compiles the source code
- Execute: mpirun -np <number of processors> ./<executable name>  <time to run (seconds)>
  
## SpeedUp
- Achieved SpeedUp ~x11 on 4 nodes - 12 processors each node. 


## Results
- Open output (.dat) file to see the results of the analytical and numerical solution.
  
## Contact
- apost035@umn.edu, trs.apostolou@gmail.com
