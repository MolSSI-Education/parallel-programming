#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
    if (MPI_Init(&argc,&argv) != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, 1);
    
    int nproc, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (nproc != 2)  MPI_Abort(MPI_COMM_WORLD, 2);

    
    int tag = 33;
    if (rank == 0) {
        int value = 99;
        MPI_Send(&value, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
    }
    else {
        int value=0;
        MPI_Status status;
        MPI_Recv(&value, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        if (value != 99) MPI_Abort(MPI_COMM_WORLD, 3);
        std::cout << "OK" << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}

