#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
    if (MPI_Init(&argc,&argv) != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, 1);
    
    int nproc, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (nproc < 2)  MPI_Abort(MPI_COMM_WORLD, 2);
    int left  = (rank-1+nproc)%nproc;
    int right = (rank+1)%nproc;
    
    int buffer=-1;
    int value = 99;
    int tag = 33;
    MPI_Status status;
    if (rank == 0) {
        MPI_Send(&value, 1, MPI_INT, right, tag, MPI_COMM_WORLD); // send to right
        MPI_Recv(&buffer, 1, MPI_INT, left, tag, MPI_COMM_WORLD, &status); // recv from left
        if (buffer != 99) MPI_Abort(MPI_COMM_WORLD, 3);
        std::cout << "OK" << std::endl;
    }
    else {
        MPI_Recv(&buffer, 1, MPI_INT, left, tag, MPI_COMM_WORLD, &status); // recv from left
        MPI_Send(&buffer, 1, MPI_INT, right, tag, MPI_COMM_WORLD); // send to right
    }
    
    MPI_Finalize();
    return 0;
}

