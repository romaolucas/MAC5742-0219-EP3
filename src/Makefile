OUTPUT=mandelbrot

IMAGE=.ppm

MPCC=mpicc
CC=gcc
CC_OPT=-std=c11

.PHONY: all
all: $(OUTPUT)_seq $(OUTPUT)_mpi

$(OUTPUT)_seq: $(OUTPUT)_seq.c
	$(CC) -o $(OUTPUT)_seq $(CC_OPT) $(OUTPUT)_seq.c

$(OUTPUT)_mpi: $(OUTPUT)_mpi.c
	$(MPCC) -o $@ $(CC_OPT) $^ 

.PHONY: clean
clean:
	rm $(OUTPUT)_seq *$(IMAGE)
