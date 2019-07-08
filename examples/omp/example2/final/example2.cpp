#include <stdio.h>
#include <omp.h>
#include <math.h>

int N = 1000000000;

int main()
{

  double start_time = omp_get_wtime();
  double start_loop;

  //declare two arrays
  int* a = new int[N];
  int* b = new int[N];

  //initialize a
  start_loop = omp_get_wtime();
#pragma omp parallel for
  for (int i=0; i<N; i++) {
    a[i] = 1.0;
  }
  printf("Initialize a time: %f\n",omp_get_wtime()-start_loop);

  //initialize b
  start_loop = omp_get_wtime();
#pragma omp parallel for
  for (int i=0; i<N; i++) {
    b[i] = 1.0 + double(i);
  }
  printf("Initialize b time: %f\n",omp_get_wtime()-start_loop);

  //add the two arrays
  start_loop = omp_get_wtime();
#pragma omp parallel for
  for (int i=0; i<N; i++) {
    a[i] = a[i] + b[i];
  }
  printf("Add arrays time: %f\n",omp_get_wtime()-start_loop);

  //average the result
  start_loop = omp_get_wtime();
  double average = 0.0;
#pragma omp parallel for reduction(+:average)
  for (int i=0; i<N; i++) {
    average += a[i];
  }
  average = average/double(N);

  printf("Average result time: %f\n",omp_get_wtime()-start_loop);
  printf("\n");

  //print the result
  printf("Average: %f\n",average);

  delete(a);
  delete(b);

  printf("Total time: %f\n",omp_get_wtime()-start_time);

  return 0;
}
