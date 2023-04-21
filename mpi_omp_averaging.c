#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>
#include "omp.h"
#include <mpi.h>

long N = 4;
int NI = 3;

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


int main(int argc, char **argv)
{   
    int P = 1, rank, tag_left = 0, tag_right = 1, tag_sum = 2;
    double *A, *A_shadow, *check, global_sum = 0.0, local_sum = 0.0;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);


    if(getenv("N") != NULL){
        N = atoi(getenv("N"));
    }
    if(getenv("NI")!=NULL){
        NI = atoi(getenv("NI"));
    }
    if(getenv("N")!=NULL){
        P = atoi(getenv("P"));
    }

    A = (double*)malloc((N+2)*sizeof(double));
    A_shadow = (double*)malloc((N+2)*sizeof(double));
    check = (double*)malloc((N+2)*sizeof(double));

    memset(A, 0, (N+2)*sizeof(double));
    memset(A_shadow, 0, (N+2)*sizeof(double));
    memset(check, 0, (N+2)*sizeof(double));

    A[N+1] = N+1;
    A_shadow[N+1] = N+1;
    MPI_Status status;

    // long chunk_size = ((N+2)/P);

    // long start = (rank*chunk_size) + 1;
    // long end = (rank+1)*chunk_size;
    // if(rank == P-1){
    //     end = N+1;
    // }
    int batchSize = ceilDiv((long)P);
    long start = rank * batchSize + 1;
    long end = min(start + batchSize - 1, N);

    long start_time = get_usecs();

    // Sending edge elements calculated by this rank
    if(rank > 0){
        MPI_Send(&A[start], 1, MPI_DOUBLE, rank-1, tag_left, MPI_COMM_WORLD);
    }
    if(rank < P-1){
        MPI_Send(&A[end], 1, MPI_DOUBLE, rank+1, tag_right, MPI_COMM_WORLD);
    }

    if(rank > 0){
        MPI_Recv(&A[start-1], 1, MPI_DOUBLE, rank-1, tag_left, MPI_COMM_WORLD, &status);
    }
    if(rank < P-1){
        MPI_Recv(&A[end+1], 1, MPI_DOUBLE, rank+1, tag_right, MPI_COMM_WORLD, &status);
    }

    // int batchSize = ceilDiv((long)P);

    for (int iter = 0; iter < NI; iter++){

        // long start = rank * batchSize + 1;
        // long end = min(start + batchSize - 1, N);
        #pragma omp parallel for 
        for (int j = start; j <= end; j++){
            A_shadow[j] = (A[j-1] + A[j+1]) / 2.0;
        }

        double* temp = A_shadow;
        A_shadow = A;
        A = temp;

    }

    // Calculating local sum of the chunk
    for(int i = start ; i <= end ; i++){
        local_sum += A[i];
    }

    // Sending local sum to master rank 0
    if(rank > 0){
        int dest = 0;
        MPI_Send(&local_sum, 1, MPI_DOUBLE, dest, tag_sum, MPI_COMM_WORLD);
    }
    // Collecting sums from all processors at rank 0 processor
    else {
        global_sum = local_sum;
        for(int i = 1; i < P ; i++) {
            int src = i;
            double recv_sum;
            MPI_Recv(&recv_sum, 1, MPI_DOUBLE, src, tag_sum, MPI_COMM_WORLD, &status);
            global_sum += recv_sum; 
        }
        printf("\nSum of array A is: %f\n", global_sum);
    }

    long end_time = get_usecs();
    double dur = ((double)(end_time - start_time))/1000000;

    printf("Time = %.5f\n",dur);
    
    free(A);
    free(A_shadow);
    free(check); 

    // Finalize MPI
    MPI_Finalize();

    return 0;
}