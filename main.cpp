#include <cstring>
#include <cstdlib>
#include <fstream>
#include "gates.h"

using namespace std;

typedef complex<double> complexd;
unsigned long long len = 0;

int handleError(int errorCode, int rank, complexd **states, complexd **statesFriend) {
	if (!rank) {
		switch (errorCode) {
			case 1: cout << "Error 01: no parameters" << endl; break;
			case 2: cout << "Error 02: invalid mode" << endl; break;
			case 3: cout << "Error 03: invalid parameters: no input file" << endl; break;
			case 4: cout << "Error 04: invalid parameters: no count of qubits" << endl; break;
			case 5: cout << "Error 05: invalid parameters: no output file" << endl; break;
			case 6: cout << "Error 06: invalid parameters: invalid input file" << endl; break;
			case 7: cout << "Error 07: invalid parameters: invalid output file" << endl; break;
			case 9: cout << "Error 09: invalid parameters: invalid write mode" << endl; break;
			case 10: cout << "Error 10: invalid parameters: no valid file" << endl; break;
			case 11: cout << "Error 11: invalid parameters: no gate" << endl; break;
			case 12: cout << "Error 12: invalid parameters: no qubit for gate" << endl; break;
			case 13: cout << "Error 13: invalid parameters: no second qubit for gate" << endl; break;
			case 14: cout << "Error 14: invalid parameters: invalid gate" << endl; break;
			case 15: cout << "Error 15: invalid parameters: invalid valid file" << endl; break;
			case 16: cout << "Error 16: invalid parameters: no time file" << endl; break;
		}
		cout << "usage:" << endl;
		cout << "\tmpirun -np <threads> go read <inputFile> <countQubits> <timeFile> <outputFile> <gate> [<qubit> <secondQubit>]" << endl;
		cout << "\tmpirun -np <threads> go generate <countQubits> <outputFile>" << endl;
		cout << "\tmpirun -np <threads> go test <inputFile> <countQubits> <validFile> <gate> [<qubit> <secondQubit>]" << endl;
	}
	MPI_Finalize();
	if (states)
		delete [] *states;
	if (statesFriend)
		delete [] *statesFriend;
	return 0;
}

int writeVectorInFile(complexd *states, char *file, int size, int rank) {
	MPI_File fd;
	MPI_Status status;
	if (MPI_File_open(MPI_COMM_WORLD, file, MPI_MODE_WRONLY, MPI_INFO_NULL, &fd) != MPI_SUCCESS)
		return 1;
	MPI_File_seek(fd, len * rank * 2 * sizeof(double), MPI_SEEK_SET);
	for (int i = 0; i < len; i++) {
		double buf[2];
		buf[0] = states[i].real();
		buf[1] = states[i].imag();
		MPI_File_write(fd, buf, 2, MPI_DOUBLE, &status);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&fd);
	return 0;
}

void generateVector(complexd *states, int size, int rank) {
	double norm = 0;
	unsigned int seedp = time(NULL) + rank;

	for (long long i = 0; i < len; i++) {
		double realPart = (double) rand_r(&seedp) / RAND_MAX;
		double imagPart = (double) rand_r(&seedp) / RAND_MAX;
		states[i] = complexd(realPart, imagPart);
		norm += realPart * realPart + imagPart * imagPart;
	}
	// нормирование вектора
	double normSum = 0;
	MPI_Allreduce(&norm, &normSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	normSum = sqrt(normSum);
	for (long long i = 0; i < len; i++)
		states[i] /= normSum;
}

int readVectorFromFile(complexd *states, char *file, int size, int rank) {
	MPI_File fd;
	MPI_Status status;
	if (MPI_File_open(MPI_COMM_WORLD, file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fd) != MPI_SUCCESS)
		return 1;
	MPI_File_seek(fd, len * rank * 2 * sizeof(double), MPI_SEEK_SET);
	for (int i = 0; i < len; i++) {
		double buf[2];
		MPI_File_read(fd, buf, 2, MPI_DOUBLE, &status);
		states[i] = complexd(buf[0], buf[1]);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&fd);
	return 0;
}

void writeVector(complexd *states, int rank) {
	for (int i = 0; i < len; i++)
		cout << states[i] << ' ' << i + rank * len << endl;
}

void testingNorm(complexd *states, int rank, int size) {
	double localSum = 0;
	for (int i = 0; i < len; i++)
		localSum += states[i].real() * states[i].real() + states[i].imag() * states[i].imag();
	double sum = 0;
	MPI_Reduce(&localSum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		if (sum < 1.000001 && sum > 0.999999)
			cout << "\t✅ blackbox" << endl;
		else
			cout << "\t❌ blackbox : " << sum << endl;
	}
}

void testing(complexd *firstStates, complexd *secondStates, int rank, int size) {
	bool result = true;
	for (int i = 0; i < len; i++)
		if (abs(firstStates[i].real() - secondStates[i].real()) > 0.000001 ||
		abs(firstStates[i].imag() - secondStates[i].imag()) > 0.000001)  {
			result = false;
			break;
		}
	bool results[size];
	MPI_Gather(&result, 1, MPI_BYTE, results, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
	if (!rank) {
		for (int i = 0; i < size; i++)
			if (results[i] == false) {
				result = false;
				break;
			}
		if (result == true)
			cout << "\t✅ canonization" << endl;
		else
			cout << "\t❌ canonization" << endl;
	}
}

int main(int argc, char **argv) {
	int size, rank;
	MPI_Status status;
	MPI_Request request;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc == 1)
		return handleError(1, rank, 0, 0);

	if (!strcmp(argv[1], "generate")) {
		switch (argc) {
			case 2: return handleError(4, rank, 0, 0);
			case 3: return handleError(5, rank, 0, 0);
		}
		int countQubits = atoi(argv[2]);
		char *outFile = argv[3];
		len = pow((double) 2, countQubits) / size;
		complexd *states = new complexd[len];
		generateVector(states, size, rank);
		if (writeVectorInFile(states, outFile, size, rank))
			return handleError(7, rank, &states, 0);
		delete [] states;
	} else if (!strcmp(argv[1], "read")) {
		switch (argc) {
			case 2: return handleError(3, rank, 0, 0);
			case 3: return handleError(4, rank, 0, 0);
			case 4: return handleError(16, rank, 0, 0);
		}
		char *inFile = argv[2];
		int countQubits = atoi(argv[3]);
		len = pow((double) 2, countQubits) / size;
		complexd *states = new complexd[len];
		generateVector(states, size, rank);
		// if (readVectorFromFile(states, inFile, size, rank))
		// 	return handleError(6, rank, &states, 0);
		char *timeFile = argv[4];
		char *outputFile = argv[5];
		char *gate = argv[6];
		double timeStart = MPI_Wtime();
		if (!strcmp(gate, "Hn"))
			Hn(&states, countQubits, rank, size);
		else {
			if (argc == 6)
				return handleError(12, rank, &states, 0);
			int qubit = atoi(argv[7]);
			if (!strcmp(gate, "H")) {
				H(&states, countQubits, qubit, rank, size);
			} else if (!strcmp(gate, "NOT")) {
				NOT(&states, countQubits, qubit, rank, size);
			} else if (!strcmp(gate, "CNOT")) {
				if (argc == 7)
					return handleError(13, rank, &states, 0);
				int secondQubit = atoi(argv[8]);
				CNOT(&states, countQubits, qubit, secondQubit, rank, size);
			} else if (!strcmp(gate, "ROT")) {
				ROT(&states, countQubits, qubit, rank, size);
			} else if (!strcmp(gate, "CROT")) {
				if (argc == 7)
					return handleError(13, rank, &states, 0);
				int secondQubit = atoi(argv[8]);
				CROT(&states, countQubits, qubit, secondQubit, rank, size);
			} else
				return handleError(14, rank, &states, 0);
		}
		double timeEnd = MPI_Wtime();
		double time = timeEnd - timeStart;
		if (strcmp(timeFile, "null")) {
			double maxTime;
			MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			if (rank == 0) {
				ofstream fout(timeFile, ios_base::app);
				fout << gate << '\t' << countQubits << '\t' << size << '\t' << maxTime << endl;
				fout.close();
			}
		}
		if (strcmp(outputFile, "null"))
			if (writeVectorInFile(states, outputFile, size, rank))
				return handleError(7, rank, &states, 0);
		delete [] states;
	} else if (!strcmp(argv[1], "test")) {
		switch (argc) {
			case 2: return handleError(3, rank, 0, 0);
			case 3: return handleError(4, rank, 0, 0);
			case 4: return handleError(10, rank, 0, 0);
			case 5: return handleError(11, rank, 0, 0);
		}
		char *inFile = argv[2];
		int countQubits = atoi(argv[3]);
		len = pow((double) 2, countQubits) / size;
		complexd *states = new complexd[len];
		if (readVectorFromFile(states, inFile, size, rank))
			return handleError(6, rank, &states, 0);
		char *validFile = argv[4];
		char *gate = argv[5];
		if (!strcmp(gate, "Hn"))
			Hn(&states, countQubits, rank, size);
		else {
			if (argc == 6)
				return handleError(12, rank, &states, 0);
			int qubit = atoi(argv[6]);
			if (!strcmp(gate, "H")) {
				H(&states, countQubits, qubit, rank, size);
			} else if (!strcmp(gate, "NOT")) {
				NOT(&states, countQubits, qubit, rank, size);
			} else if (!strcmp(gate, "CNOT")) {
				if (argc == 7)
					return handleError(13, rank, &states, 0);
				int secondQubit = atoi(argv[7]);
				CNOT(&states, countQubits, qubit, secondQubit, rank, size);
			} else if (!strcmp(gate, "ROT")) {
				ROT(&states, countQubits, qubit, rank, size);
			} else if (!strcmp(gate, "CROT")) {
				if (argc == 7)
					return handleError(13, rank, &states, 0);
				int secondQubit = atoi(argv[7]);
				CROT(&states, countQubits, qubit, secondQubit, rank, size);
			} else
				return handleError(14, rank, &states, 0);
		}
		testingNorm(states, rank, size);
		complexd *validStates = new complexd[len];
		if (readVectorFromFile(validStates, validFile, size, rank))
			return handleError(15, rank, &states, &validStates);
		testing(states, validStates, rank, size);
		delete [] validStates;
		delete [] states;

	} else
		return handleError(2, rank, 0, 0);

	MPI_Finalize();
	return 0;
}
