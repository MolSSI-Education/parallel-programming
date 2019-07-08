#include <mpi.h>
#include <cstdio>
using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    const int maxlen = 1024*1024; // Vary this ... 1, 2, 128, 1024, 1024*1024
    int nproc, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int other = (rank+1)%2;

    if (nproc != 2) MPI_Abort(MPI_COMM_WORLD, 1);

    char* buf1 = new char[maxlen];
    char* buf2 = new char[maxlen];

    // This is the unsafe approach
    // if (rank == 0) printf("trying 1\n");
    // MPI_Send(buf1, maxlen, MPI::BYTE, other, 1, MPI_COMM_WORLD);
    // MPI_Recv(buf2, maxlen, MPI::BYTE, other, 1, MPI_COMM_WORLD, 
    // 	     MPI_STATUS_IGNORE);

    if (rank == 0) printf("trying 2\n");
    if (rank == 0) {
      MPI_Send(buf1, maxlen, MPI::BYTE, other, 1, MPI_COMM_WORLD);
      MPI_Recv(buf2, maxlen, MPI::BYTE, other, 1, MPI_COMM_WORLD, 
	       MPI_STATUS_IGNORE);
    }
    else {
      MPI_Recv(buf2, maxlen, MPI::BYTE, other, 1, MPI_COMM_WORLD, 
	       MPI_STATUS_IGNORE);
      MPI_Send(buf1, maxlen, MPI::BYTE, other, 1, MPI_COMM_WORLD);
    }

    if (rank == 0) printf("trying 3\n");
    MPI_Request req;
    MPI_Irecv(buf2, maxlen, MPI::BYTE, other, 3, MPI_COMM_WORLD, &req);
    MPI_Send(buf1, maxlen, MPI::BYTE, other, 3, MPI_COMM_WORLD);
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    if (rank == 0) printf("trying 3a\n");
    MPI_Isend(buf1, maxlen, MPI::BYTE, other, 3, MPI_COMM_WORLD, &req);
    MPI_Recv(buf2, maxlen, MPI::BYTE, other, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Wait(&req, MPI_STATUS_IGNORE);

    if (rank == 0) printf("trying 4\n");
    MPI_Sendrecv(buf1, maxlen, MPI::BYTE, other, 2,
		 buf2, maxlen, MPI::BYTE, other, 2, 
		 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    delete buf1;
    delete buf2;

    MPI_Finalize();
    return 0;
}
