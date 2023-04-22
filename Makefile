EXE = omp_averaging sequential_averaging mpi_omp_averaging mm_mpi

all: $(EXE)

omp_averaging: omp_averaging.c
	gcc -O3 -o $@ $^ -fopenmp

mpi_omp_averaging: mpi_omp_averaging.c
	mpicc -O3 -o $@ $^ -fopenmp

mm_mpi: mm_mpi.c
	mpicc -O3 -o $@ $^ 

sequential_averaging: sequential_averaging.c
	gcc -O3 -o $@ $^

clean:
	rm -rf *.o $(EXE)