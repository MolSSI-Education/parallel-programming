#include <mpi.h>
#include <iostream>
#include <cmath>
#include <cstdlib>

double drand() {
    const double fac = 1.0/(RAND_MAX-1.0);
    return fac*random();
}

void set_seed(int p) {
    srandom(p);
    for (int i=0; i<100; i++) drand();
}

int main(int argc, char** argv) {
    MPI_Init(&argc,&argv);
    int nproc, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    set_seed(rank);
    
    const long N = 100000000/nproc;
    long sum = 0;
    for (long i=0; i<N; i++) {
        double x = 2.0*(drand()-0.5); // Random value in [-1,1]
        double y = 2.0*(drand()-0.5); // Random value in [-1,1]
        double rsq = x*x + y*y;
        if (rsq < 1.0) sum++;
    }
    long total;
    MPI_Reduce(&sum, &total, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double pi = (4.0*total)/(N*nproc);
        std::cout.precision(8);
        std::cout << pi << std::endl;
    }

    MPI_Finalize();
    return 0;
}

