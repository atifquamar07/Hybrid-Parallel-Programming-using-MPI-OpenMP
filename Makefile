EXE = omp_averaging sequential_averaging mpi_omp_averaging mpi_omp_averaging_blocking mm_mpi mm_mpi_blocking

all: $(EXE)

omp_averaging: omp_averaging.c
	gcc -o $@ $^ -fopenmp

mpi_omp_averaging: mpi_omp_averaging.c
	mpicc -o $@ $^ -fopenmp

mpi_omp_averaging_blocking: mpi_omp_averaging_blocking.c
	mpicc -o $@ $^ -fopenmp

mm_mpi: mm_mpi.c
	mpicc -o $@ $^ 

mm_mpi_blocking: mm_mpi_blocking.c
	mpicc -o $@ $^ 

sequential_averaging: sequential_averaging.c
	gcc -o $@ $^

clean:
	rm -rf *.o $(EXE)