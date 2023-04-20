#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>
#include "omp.h"
#include <mpi.h>


long N = 51658240;
int NI = 30, P = 1, rank, tag = 42;
double *A, *A_shadow, *check, global_sum = 0.0;

long get_usecs () {

    struct timeval t;
    gettimeofday(&t,NULL);
    return t.tv_sec*1000000 + t.tv_usec;

}

void runParallel(){

    for (int iter = 0; iter < NI; iter++){

        int chunk_size = ((N+2)/P);

        int start = (rank*chunk_size) + 1;
        int end = (rank+1)*chunk_size;
        if(rank == P-1){
            end = N+1;
        }


        MPI_Status status;

        if(rank > 0){
            MPI_Recv(&A[start-1], 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
        }
        if(rank < P-1){
            MPI_Recv(&A[end+1], 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
        }

        // #pragma omp parallel for default(none) firstprivate(start, end) shared(A, A_shadow, rank, tag) 
        for (int j = start; j <= end; j++) {
            A_shadow[j] = (A[j-1] + A[j+1]) / 2.0;
            if(j == start && rank > 0){
                MPI_Send(&A_shadow[start], 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
            }
        }
        if(rank < P-1){
            MPI_Send(&A_shadow[end], 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
        }
        
        double* temp = A_shadow;
        A_shadow = A;
        A = temp;

    }

}

// double arraySum(){

//     double local_sum = 0.0;

//     int start = (rank*chunk_size) + 1;
//     int end = (rank+1)*chunk;
//     if(end > N+1){
//         end = N+1;
//     }

//     #pragma omp parallel for default(none) firstprivate(start, end) shared(A) reduction(+:sum)
//     for(int i = start ; i <= end ; i++){
//         sum += A[i];
//     }

//     return sum;
// }


int main(int argc, char **argv)
{   
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

    long start = get_usecs();
    runParallel();
    long end = get_usecs();
    double dur = ((double)(end-start))/1000000;

    printf("Time = %.5f\n",dur);
    
    // start = get_usecs();
    // double sum = arraySum();
    // end = get_usecs();
    // dur = ((double)(end-start))/1000000;

    // printf("\nSum of array A is: %f\n", sum);
    // printf("Time = %.5f\n",dur);

    free(A);
    free(A_shadow);
    free(check); 

    // Finalize MPI
    MPI_Finalize();


    return 0;
}