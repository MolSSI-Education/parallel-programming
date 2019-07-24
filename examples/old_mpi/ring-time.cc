#include <mpi.h>
#include <cstdio>
#include <algorithm>
using namespace std;

int main(int argc, char** argv) {
    const int MAXLEN = 1024*1024;
    char* buf = new char[MAXLEN];

    int nproc, rank;
    MPI::Init(argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int left  = (rank+nproc-1)%nproc;
    const int right = (rank+1)%nproc;

    if (nproc < 2) MPI_Abort(MPI_COMM_WORLD, 1);

    if (rank == 0) {
        printf("   msglen     nloop      used    rate (byte/s)\n");
        printf("  -------  --------  --------  --------\n");
    }

    for (int msglen=1; msglen<=MAXLEN; msglen*=2) {
        double testim = (5e-6 + msglen*1e-9)*nproc;
        int nloop = max(0.1/testim,1.0);

        MPI_Barrier(MPI_COMM_WORLD);
        double used = MPI::Wtime();
        for (int loop=0; loop<nloop; loop++) {
            if (rank == 0) {
	      MPI_Send(buf, msglen, MPI::BYTE, right, 1, MPI_COMM_WORLD);
	      MPI_Recv(buf, msglen, MPI::BYTE, left,  1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else {
	      MPI_Recv(buf, msglen, MPI::BYTE, left,  1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      MPI_Send(buf, msglen, MPI::BYTE, right, 1, MPI_COMM_WORLD);
            }
        }
        used = MPI::Wtime() - used;
        double rate = nloop*nproc*msglen/used;
        if (rank == 0) printf(" %8d  %8d  %.2e   %.2e\n", msglen, nloop, used, rate);
    }

    delete buf;

    MPI::Finalize();
    return 0;
}
