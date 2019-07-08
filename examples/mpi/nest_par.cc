#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc,&argv);
    
    int nproc, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int count = 0;
    double sum = 0.0;
    for (int i=0; i<10; i++) {
      for (int j=0; j<i; j++, count++) {
	if ((count%nproc)==rank) {
	  sum += j;
	}
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) std::cout << sum << std::endl;

    MPI_Finalize();
    return 0;
}

