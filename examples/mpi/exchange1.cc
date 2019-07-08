#include <mpi.h>
#include <iostream>
#include <cstdlib>

int main(int argc, char** argv) {
    if (MPI_Init(&argc,&argv) != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, 1);
    int rank, nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (nproc != 2)  MPI_Abort(MPI_COMM_WORLD, 2);

    size_t N=0;
    if (argc == 2) N = atoi(argv[1]);
    if (argc != 2 || N < 1 || N>100000000) {
        std::cerr << "usage: 'exchange1 N' where N is a positive integer < 10**8" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
    if (rank == 0) std::cout << "N " << N << std::endl;
    
    int other = (rank+1)%nproc;
    int tag = 3;

    unsigned char* buf0 = new unsigned char[N];
    unsigned char* buf1 = new unsigned char[N];

    for (size_t i=0; i<N; i++) buf0[i] = (unsigned char)(rank);
    for (size_t i=0; i<N; i++) buf1[i] = 99;

    MPI_Status status;
    MPI_Send(buf0, N, MPI_UNSIGNED_CHAR, other, tag, MPI_COMM_WORLD);
    MPI_Recv(buf1, N, MPI_UNSIGNED_CHAR, other, tag, MPI_COMM_WORLD, &status);
    
    for (size_t i=0; i<N; i++) {
        if (buf1[i] != other) {
            std::cerr << "bad value " << i << " " << buf1[i] << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 3);
        }
    }

    if (rank == 0) std::cout << "OK" << std::endl;
    
    MPI_Finalize();
    return 0;
}

