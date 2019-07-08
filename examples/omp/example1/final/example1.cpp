#include <stdio.h>
#include <omp.h>

int main()
{
#pragma omp parallel
  {
    int thread_id = omp_get_thread_num();
    printf("Hello World! (%i)\n", thread_id);
  }

  int thread_id = omp_get_thread_num();
  printf("This is another message! (%i)\n", thread_id);

#pragma omp parallel num_threads(3)
  {
    int thread_id = omp_get_thread_num();
    printf("Goodbye World! (%i)\n", thread_id);
  }

  return 0;
}
