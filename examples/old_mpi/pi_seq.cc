#include <iostream>
#include <cmath>
#include <cstdlib>

double drand() {
    const double fac = 1.0/(RAND_MAX-1.0);
    return fac*random();
}

int main() {
  const long N = 100000000;
  long sum = 0;
  for (long i=0; i<N; i++) {
    double x = 2.0*(drand()-0.5); // Random value in [-1,1]
    double y = 2.0*(drand()-0.5); // Random value in [-1,1]
    double rsq = x*x + y*y;
    if (rsq < 1.0) sum++;
  }
  double pi = (4.0*sum)/N;
  std::cout.precision(8);
  std::cout << pi << std::endl;
  return 0;
}

