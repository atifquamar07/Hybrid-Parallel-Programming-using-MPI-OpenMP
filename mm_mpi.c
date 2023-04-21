#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>
// #include <mpi.h>

int N = 500;

long get_usecs () {
    struct timeval t;
    gettimeofday(&t,NULL);
    return t.tv_sec*1000000 + t.tv_usec;
}

long ceilDiv(long d) {
    long m = N / d;
    if (m * d == N) {
        return m;
    } else {
        return (m + 1);
    }
}

long min(long a, long b){
    if(a > b){
        return b;
    }
    else {
        return a;
    }
}


int main(int argc, char **argv){

    int P = 1, rank, tag = 10;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    int **A = (int **)malloc(N * sizeof(int *)); 
    int **B = (int **)malloc(N * sizeof(int *)); 
    int **C = (int **)malloc(N * sizeof(int *)); 

    for (int i = 0; i < N; i++) {
        A[i] = (int *)malloc(N * sizeof(int)); 
        B[i] = (int *)malloc(N * sizeof(int)); 
        C[i] = (int *)malloc(N * sizeof(int));   
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = 1;
        }
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            B[i][j] = 1;
        }
    }

    long start_time = get_usecs();

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    long end_time = get_usecs();
    double dur = ((double)(end_time - start_time))/1000000;

    printf("Time = %.5f\n",dur);

}