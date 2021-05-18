G++ = mpic++
FLAG = -O2 -std=c++11
NAME = go
SRCS = main.cpp

.PHONY: all test generate fclean report

all: compile
	@echo "usage:"
	@echo "	mpirun -np <threads> go read <inputFile> <countQubits> <timeFile> <outputFile> <gate> [<qubit> <secondQubit>]"
	@echo "	mpirun -np <threads> go generate <countQubits> <outputFile>"
	@echo "	mpirun -np <threads> go test <inputFile> <countQubits> <validFile> <gate> [<qubit> <secondQubit>]"

compile:
	@$(G++) $(FLAG) $(SRCS) -o $(NAME)
	@echo "Compiled"

time: compile
	for qubits in 25 26 27 ; do \
		for threads in 1 2 4 8 16 32 64 ; do \
			mpisubmit.pl -p $$threads go read null $$qubits time.txt null Hn ; \
			mpisubmit.pl -p $$threads go read null $$qubits time.txt null CNOT 1 2 ; \
	done ; done

generate: compile
	@touch test/vector.txt
	@echo "Generating vector (16 qubits)"
	@mpirun -np 8 go generate 16 test/vector.txt
	@for GATE in H NOT ROT ; do \
		for num in 1 2 16 ; do \
			touch "test/$$GATE $$num.txt" ; \
			echo "generating gate $$GATE with qubit $$num" ; \
			mpirun -np 8 go read test/vector.txt 16 null "test/$$GATE $$num.txt" $$GATE $$num ; \
	done ; done
	@touch "test/Hn.txt"
	@echo "generating gate Hn"
	@mpirun -np 8 go read test/vector.txt 16 null "test/Hn.txt" Hn
	@for GATE in CNOT CROT ; do \
		for nums in "1 2" "2 8" "15 16" ; do \
			touch "test/$$GATE $$nums.txt" ; \
			echo "generating gate $$GATE with qubits $$nums" ; \
			mpirun -np 8 go read test/vector.txt 16 null "test/$$GATE $$nums.txt" $$GATE $$nums ; \
	done ; done


test: compile
	@for GATE in H NOT ROT ; do \
		for num in 1 2 16 ; do \
			echo "testing gate $$GATE with qubit $$num" ; \
			mpirun -np 8 go test test/vector.txt 16 "test/$$GATE $$num.txt" $$GATE $$num ; \
	done ; done
	@echo "testing gate Hn"
	@mpirun -np 8 go test test/vector.txt 16 "test/Hn.txt" Hn
	@for GATE in CNOT CROT ; do \
		for nums in "1 2" "2 8" "15 16" ; do \
			echo "testing gate $$GATE with qubits $$nums" ; \
			mpirun -np 8 go test test/vector.txt 16 "test/$$GATE $$nums.txt" $$GATE $$nums ; \
	done ; done

module:
	@echo "module load SpectrumMPI/10.1.0"

clean:
	@rm $(NAME)
	@echo "Cleaned"
