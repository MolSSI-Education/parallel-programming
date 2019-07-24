#include <cstdio>
#include <cmath>

double g(double x) {
  return exp(-x*x)*cos(3*x);
}


// Make this routine run in parallel
double integrate(const int N, const double a, const double b, double (*f)(double)) {
  const double h = (b-a)/N;
  double sum=0.0;
  for (int i=1; i<(N-1); i++) {
    sum += f(a + i*h);
  }
  sum += 0.5*(f(b) + f(a));
  return sum*h;
}

int main(int argc, char** argv) {
  const double gexact = std::sqrt(4.0*std::atan(1.0))*std::exp(-9.0/4.0);
  const double a=-6.0, b=6.0;

  double result_prev = integrate(1, a, b, g);

  // Process 0 should communicate N to all other processes
  for (int N=2; N<=1024; N*=2) {
    double result = integrate(N, a, b, g);
    double err_prev = std::abs(result-result_prev);
    double err_exact = std::abs(result-gexact);
    printf("N=%2d   result=%.10e   err-prev=%.2e   err-exact=%.2e\n",
	   N, result, err_prev, err_exact);

    // Please have only process 0 determine if we are converged and somehow
    // tell everyone else we are finished.
    if (err_prev < 1e-10 && N>4) break;
    result_prev = result;
  }
   
  return 0;
}
