#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>
#include <mpi.h>

int N = 1024;

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

    int P = 1, rank, tag_a = 1, tag_b = 2, tag_start = 3, tag_end = 4, tag_index = 3, tag_rows = 4;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Status status;

    int **A = (int **)malloc(N * sizeof(int *)); 
    int **B = (int **)malloc(N * sizeof(int *)); 
    int **C = (int **)malloc(N * sizeof(int *)); 

    for (int i = 0; i < N; i++) {
        A[i] = (int *)malloc(N * sizeof(int)); 
        B[i] = (int *)malloc(N * sizeof(int)); 
        C[i] = (int *)malloc(N * sizeof(int));   
    }

    
    if(rank == 0){
        
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

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                C[i][j] = -1;
            }
        }

        long start_time = get_usecs();

        if(P > 1){

            for (int dest = 1; dest < P; dest++) {
                for(int i = 0 ; i < N ; i++){
                    MPI_Send(&A[i][0], N, MPI_INT, dest, tag_a, MPI_COMM_WORLD);
                    MPI_Send(&B[i][0], N, MPI_INT, dest, tag_b, MPI_COMM_WORLD);
                }
            }
            // int source = 0;
           
            // for(int i = 0 ; i < N ; i++){
            //     MPI_Bcast(&A[i][0], N, MPI_INT, source, MPI_COMM_WORLD);
            //     MPI_Bcast(&B[i][0], N, MPI_INT, source, MPI_COMM_WORLD);
            // }
            

            for(int i = 1 ; i < P ; i++){
                int start_master, end_master;
                MPI_Recv(&start_master, 1, MPI_INT, i, tag_start, MPI_COMM_WORLD, &status);
                MPI_Recv(&end_master, 1, MPI_INT, i, tag_end, MPI_COMM_WORLD, &status);

                for (int j = start_master; j <= end_master; j++) {
                    MPI_Recv(&C[j][0], N, MPI_INT, i, tag_rows, MPI_COMM_WORLD, &status);
                }
            }

        }
        

        int batchSize = ceilDiv((long)P);
        int start = (rank * batchSize);
        int end = min(start + batchSize - 1, N-1);

        for (int i = start; i <= end; ++i) {
            for (int j = 0; j < N; ++j) {
                C[i][j] = 0;
                for (int k = 0; k < N; ++k) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }


        long end_time = get_usecs();
        double dur = ((double)(end_time - start_time))/1000000;
        



        // Validating output
        int **CHECK = (int **)malloc(N * sizeof(int *)); 
        for (int i = 0; i < N; i++){
            CHECK[i] = (int *)malloc(N * sizeof(int)); 
        }

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                CHECK[i][j] = 0;
                for (int k = 0; k < N; ++k) {
                    CHECK[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        int ok = 1;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if(C[i][j] != CHECK[i][j]){
                    printf("ERROR: Validation Failed!\n");
                    printf("Difference: CHECK[%d][%d] = %d != C[%d][%d] = %d\n", i, j, CHECK[i][j], i, j, C[i][j]);
                    ok = 0;
                }
            }
        }
        if(ok == 1){
            printf("OK.\n");
            printf("N = %d\n", N);
            printf("Processor Count: %d\n", P);
            printf("Time = %.5f\n",dur);
        }

    }
    if(rank > 0) {

        int from = 0;

        for(int i = 0 ; i < N ; i++){
            MPI_Recv(&A[i][0], N, MPI_INT, from, tag_a, MPI_COMM_WORLD, &status);
            MPI_Recv(&B[i][0], N, MPI_INT, from, tag_b, MPI_COMM_WORLD, &status);
        }
   
        
        // printf("In rank %d, matrices A and B received are: \n", rank);
        // for (int i = 0; i < N; ++i) {
        //     for (int j = 0; j < N; ++j) {
        //         printf("%d ", A[i][j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");

        // for (int i = 0; i < N; ++i) {
        //     for (int j = 0; j < N; ++j) {
        //         printf("%d ", A[i][j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");

        int batchSize = ceilDiv((long)P);
        int start = (rank * batchSize) ;
        int end = min(start + batchSize - 1, N-1);

        for (int i = start; i <= end; ++i) {
            for (int j = 0; j < N; ++j) {
                C[i][j] = 0;
                for (int k = 0; k < N; ++k) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        int dest = 0;
        MPI_Send(&start, 1, MPI_INT, dest, tag_start, MPI_COMM_WORLD);
        MPI_Send(&end, 1, MPI_INT, dest, tag_end, MPI_COMM_WORLD);

        for (int i = start; i <= end; i++) {
            MPI_Send(&C[i][0], N, MPI_INT, dest, tag_rows, MPI_COMM_WORLD);
        }

    }


    free(A);
    free(B);
    free(C); 
    

    // Finalize MPI
    MPI_Finalize();

}