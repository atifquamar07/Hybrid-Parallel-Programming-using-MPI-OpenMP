#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>
#include "omp.h"
#include <mpi.h>

long N = 25165824;
int NI = 64;

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
    long start_time, end_time;
  
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
    if(getenv("P")!=NULL){
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

    int batchSize = ceilDiv((long)P);
    long start = (rank * batchSize) + 1;
    long end = min(start + batchSize - 1, N);

    if(rank == 0){
        start_time = get_usecs();
    }
    

    // Sending edge elements calculated by this rank
    if(rank > 0){
        MPI_Send(&A[start], 1, MPI_DOUBLE, (rank-1), tag_left, MPI_COMM_WORLD);
    }
    if(rank < (P-1)){
        MPI_Send(&A[end], 1, MPI_DOUBLE, (rank+1), tag_right, MPI_COMM_WORLD);
    }

    if(rank > 0){
        MPI_Recv(&A[start-1], 1, MPI_DOUBLE, (rank-1), tag_right, MPI_COMM_WORLD, &status);
    }
    if(rank < (P-1)){
        MPI_Recv(&A[end+1], 1, MPI_DOUBLE, (rank+1), tag_left, MPI_COMM_WORLD, &status);
    }

    for (int iter = 0; iter < NI; iter++){

        #pragma omp parallel for default(none) firstprivate(start, end) shared(A, A_shadow) 
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

    // COllecting sums from all processors using MPI_Reduce
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0){
        end_time = get_usecs();
        double dur = ((double)(end_time - start_time))/1000000;
        printf("Time = %.5f\n",dur);
        printf("\nSum of array A is: %f\n", global_sum);
    }
    
    
    free(A);
    free(A_shadow);
    free(check); 

    // Finalize MPI
    MPI_Finalize();

    return 0;
}