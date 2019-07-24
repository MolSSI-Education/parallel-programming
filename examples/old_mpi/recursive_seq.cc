#include <cstdio>
#include <cmath>
#include <algorithm>

double g(double x) {
  return exp(-x*x)*cos(3*x);
}

// 3-point Gauss-Legendre quadrature on [a,b]
double I(const double a, const double b, double (*f)(double)) {
  // Points/weights for GL on [0,1]
  const int N = 3;
  const double x[3] = {8.87298334620741702e-01, 5.00000000000000000e-01, 1.12701665379258312e-01};
  const double w[3] = {2.77777777777777790e-01, 4.44444444444444420e-01, 2.77777777777777790e-01};

  double L = (b-a);
  double sum = 0.0;
  for (int i=0; i<N; i++) {
    sum += w[i]*f(a+L*x[i]);
  }
  return sum*L;
}

double Irecur(const double a, const double b, double (*f)(double), const double eps, int level=0) {
  const double middle = (a+b)*0.5;
  double total=I(a,b,f);
  double left = I(a,middle,f);
  double right= I(middle,b,f);
  double test=left+right;
  double err=std::abs(total-test);

  for (int i=0; i<level; i++) printf("  ");
  printf("[%6.3f,%6.3f] total=%.6e test=%.6e err=%.2e\n", a, b, total, test, err);
  if (level >= 20) return test; // 2^20 = 1M boxes

  if (err<eps) 
    return test;
  else {
    double neweps = std::max(eps*0.5,1e-15*std::abs(total));
    return Irecur(a,middle,f,neweps,level+1) + Irecur(middle,b,f,neweps,level+1);
  }
}

int main(int argc, char** argv) {
  const double gexact = std::sqrt(4.0*std::atan(1.0))*std::exp(-9.0/4.0);
  const double a=-6.0, b=6.0;
  const double eps=1e-10;

  double result = Irecur(a,b,g,eps);
  double err_exact = std::abs(result-gexact);
    printf("result=%.10e   err-exact=%.2e\n",
	   result, err_exact);

  return 0;
}
