#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <stdlib.h>
#include<sys/time.h>

// #define SIZE 25165824
// #define ITERATIONS 64

long size = 25165824;
int iter = 64;

double* A_shadow, *A;

long get_usecs () {
  struct timeval t;
  gettimeofday(&t,NULL);
  return t.tv_sec*1000000+t.tv_usec;
}

void runParallel() {
  for (int i = 0; i < iter; i++) {
    #pragma omp parallel for default(none) firstprivate(size) shared(A, A_shadow) 
    for(int j=1; j<=size; j++) {
        A_shadow[j] = (A[j - 1] + A[j + 1]) / 2.0;
    }
    double* temp = A_shadow;
    A_shadow = A;
    A = temp;
  }
}

int main(int argc, char const *argv[])
{   
    size =  25165824;
    iter = 64;
    A_shadow = (double*) malloc(sizeof(double) * (size + 2));
    A = (double*) malloc(sizeof(double) * (size + 2));

    memset(A_shadow, 0, sizeof(double) * (size + 2));
    memset(A, 0, sizeof(double) * (size + 2));
    A[size + 1] = size + 1;

    long start_time = get_usecs();
    runParallel();
    long end_time = get_usecs();
    double dur = ((double)(end_time-start_time))/1000000;
    printf("Time = %.3f\n",dur);

    double sum = 0.0;
    for(int i=1; i<=size; ++i) {
        sum+=A[i];
    }

    printf("SUM = %f\n", sum);

    free(A);
    free(A_shadow);

    return 0;
}