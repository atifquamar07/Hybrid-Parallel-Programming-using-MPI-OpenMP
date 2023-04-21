EXE = omp_averaging sequential_averaging mpi_omp_averaging

all: $(EXE)

omp_averaging: omp_averaging.c
	gcc -O3 -o $@ $^ -fopenmp

mpi_omp_averaging: mpi_omp_averaging.c
	mpicc -O3 -o $@ $^ -fopenmp

sequential_averaging: sequential_averaging.c
	gcc -O3 -o $@ $^

clean:
	rm -rf *.o $(EXE)