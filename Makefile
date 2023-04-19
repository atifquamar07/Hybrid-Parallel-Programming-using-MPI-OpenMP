EXE = omp_averaging sequential_averaging

all: $(EXE)

omp_averaging: omp_averaging.c
	gcc -O3 -o $@ $^ -fopenmp

sequential_averaging: sequential_averaging.c
	gcc -O3 -o $@ $^

clean:
	rm -rf *.o $(EXE)