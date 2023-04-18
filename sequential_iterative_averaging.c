#include<stdio.h>
#include<string.h>
#include<sys/time.h>
#include <stdlib.h>

long N = 1000000;
int NI = 10, P = 1;
double *A, *A_shadow, *myVal;

long get_usecs () {

    struct timeval t;
    gettimeofday(&t,NULL);
    return t.tv_sec*1000000 + t.tv_usec;

}

void runSequential(){

    for (int iter = 0; iter < NI; iter++){
        for (int j = 1; j <= N; j++){
            A_shadow[j] = (A[j-1] + A[j+1]) / 2.0;
        }
        double* temp = A_shadow;
        A_shadow = A;
        A = temp;
    }

}

void validateOutput() {

    int ok = 1;
    for (int i = 0; i < N + 2; i++) {
        double init = A_shadow[i];
        double curr = myVal[i];
        double diff = abs(init - curr);
        if (diff > 1e-20) {
            printf("ERROR: validation failed!\n");
            printf("Diff: myVal[%d]=%.3f != A_shadow[%d]=%.3f\n",i,curr,i,init);
            ok = 0;
            break;
        }
    }
    if(ok) printf("OK\n");

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
    myVal = (double*)malloc((N+2)*sizeof(double));

    memset(A, 0, (N+2)*sizeof(double));
    memset(A_shadow, 0, (N+2)*sizeof(double));
    memset(myVal, 0, (N+2)*sizeof(double));

    A[0] = 1.0;
    A[N+1] = 1.0;

    long start = get_usecs();
    runSequential();
    long end = get_usecs();
    double dur = ((double)(end-start))/1000000;

    printf("Time = %.3f\n",dur);
    validateOutput();
    free(A);
    free(A_shadow);
    free(myVal); 


    return 0;
}


