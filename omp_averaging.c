#include<stdio.h>
#include<string.h>
#include<sys/time.h>
#include <stdlib.h>
#include "omp.h"


long N = 25165824;
int NI = 64, P = 1, rank;
double *A, *A_shadow, *check;

long get_usecs () {

    struct timeval t;
    gettimeofday(&t,NULL);
    return t.tv_sec*1000000 + t.tv_usec;

}

void runParallel(){

    for (int iter = 0; iter < NI; iter++){

        #pragma omp parallel for default(none) firstprivate(N) shared(A, A_shadow) 
        for (int j = 1; j <= N; j++){
            A_shadow[j] = (A[j-1] + A[j+1]) / 2.0;
        }

        double* temp = A_shadow;
        A_shadow = A;
        A = temp;

    }

}

double arraySum(){

    double sum = 0.0;

    #pragma omp parallel for default(none) shared(A, N) reduction(+:sum)
    for(int i = 0 ; i < N+2 ; i++){
        sum += A[i];
    }

    return sum;
}


int main(int argc, char const *argv[])
{
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
    
    start = get_usecs();
    double sum = arraySum();
    end = get_usecs();
    dur = ((double)(end-start))/1000000;

    printf("\nSum of array A is: %f\n", sum);
    printf("Time = %.5f\n",dur);

    free(A);
    free(A_shadow);
    free(check); 


    return 0;
}