G++ = mpixlcxx
NAME = go
SRCS = main.cxx

.PHONY: all test generate fclean report

all: compile
	@echo "usage:"
	@echo "    mpirun -np <threads> go read <inputFile> <countQubits> <numQubit> <outputFile>"
	@echo "    mpirun -np <threads> go generate <outputFile> <countQubits>"
	@echo "    mpirun -np <threads> go test <inputFile> <countQubits> <numQubit> <validFile>"

compile:
	$(G++) -O2 -qstrict $(SRCS) -o $(NAME)
	@echo "Compiled"

time: compile
	touch time.txt
	for threads in 8 16 32 64 128 256 ; do \
		for j in 1 2 3 4 ; do \
			mpisubmit.bg -n $$threads go read null 28 0.01 null time.txt 1 ; \
	done ; done

test: compile
	for ((i = 0; i < 10; i++)) ; do \
		mpisubmit.bg -n 64 go test null 26 0.001 null null prec_26_0_001.txt ; \
	done ;

module:
	@echo "module load SpectrumMPI/10.1.0"

clean:
	@rm $(NAME)
	@echo "Cleaned"
