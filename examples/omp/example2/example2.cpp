#include <stdio.h>
#include <omp.h>
#include <math.h>

int N = 1000000000;

int main()
{

  //declare two arrays
  int* a = new int[N];
  int* b = new int[N];

  //initialize a
  for (int i=0; i<N; i++) {
    a[i] = 1.0;
  }

  //initialize b
  for (int i=0; i<N; i++) {
    b[i] = 1.0 + double(i);
  }

  //add the two arrays
  for (int i=0; i<N; i++) {
    a[i] = a[i] + b[i];
  }

  //average the result
  double average = 0.0;
  for (int i=0; i<N; i++) {
    average += a[i];
  }
  average = average/double(N);

  //print the result
  printf("Average: %f\n",average);

  delete(a);
  delete(b);

  return 0;
}
