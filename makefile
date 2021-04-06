G++ = mpic++
FLAG = -std=c++11
NAME = go
SRCS = main.cpp

.PHONY: all test generate fclean report

all: compile
	@echo "usage:"
	@echo "    mpirun -np <threads> go read <inputFile> <countQubits> <numQubit> <outputFile>"
	@echo "    mpirun -np <threads> go generate <outputFile> <countQubits>"
	@echo "    mpirun -np <threads> go test <inputFile> <countQubits> <numQubit> <validFile>"

compile:
	@$(G++) $(FLAG) $(SRCS) -o $(NAME)
	@echo "Compiled"

time:
	for count in 25 26 27 ; do \
	for threads in 1 2 4 8 16 32 64 ; do \
	for num in 1 2 $$count ; do \
	mpisubmit.pl -p $$threads go read vector.txt $$count $$num time.txt ; \
	done ; done ; done  

test: compile
	@touch vector.txt
	@echo "Generating vector (16 qubits) on 8 threads..."
	@mpirun -np 8 go generate vector.txt 16
	@touch 1.txt
	@for num in 1 2 16 ; do \
	echo "Testing number of qubit = $$num" ; \
	echo "Transforming vector on 1 thread..." ; \
	mpirun -np 1 go read vector.txt 16 $$num 1.txt ; \
	echo "Testing:" ; \
	for threads in 2 4 8 ; do \
	mpirun -np $$threads go test vector.txt 16 $$num 1.txt ; \
	done ; done
	@rm vector.txt 1.txt

module:
	@echo "module load SpectrumMPI/10.1.0"

clean:
	@rm $(NAME)
	@echo "Cleaned"