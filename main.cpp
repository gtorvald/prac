#include <iostream>
#include <vector>
#include <complex>
#include <cstring>
#include <cstdlib>
#include <random>
#include <fstream>
#include <cmath>
#include <mpi.h>

using namespace std;

typedef complex<double> complexd;
unsigned long long len = 0;

int handleError(int errorCode, int rank, complexd **states, complexd **statesFriend) {
	if (!rank) {
		switch (errorCode) {
			case 1: cout << "Error 01: no parameters" << endl; break;
			case 2: cout << "Error 02: invalid mode" << endl; break;
			case 3: cout << "Error 03: invalid parameters: no inputFile" << endl; break;
			case 4: cout << "Error 04: invalid parameters: no count of qubits" << endl; break;
			case 5: cout << "Error 05: invalid inputFile" << endl; break;
			case 6: cout << "Error 06: invalid parameters: no number of qubit" << endl; break;
			case 7: cout << "Error 07: invalid parameters: no outputFile" << endl; break;
			case 8: cout << "Error 08: invalid outputFile" << endl; break;
			case 9: cout << "Error 09: invalid parameters: no file with valid states" << endl; break;
			case 10: cout << "Error 10: invalid file with valid states" << endl; break;
		}
		cout << "usage:" << endl;
		cout << "\tmpirun -np <threads> go read <inputFile> <countQubits> <numQubit> <outputFile>" << endl;
		cout << "\tmpirun -np <threads> go generate <outputFile> <countQubits>" << endl;
		cout << "\tmpirun -np <threads> go test <inputFile> <countQubits> <numQubit> <validFile>" << endl;
	}
	MPI_Finalize();
	if (states)
		delete [] *states;
	if (statesFriend)
		delete [] *statesFriend;
	return 0;
}

void generateVector(complexd *states, int size, int rank) {
	mt19937 mersenne(random_device().operator()()); // генератор случайных чисел на основе Вихря Мерсенна, C++11
	double norm = 0;

	for (long long i = 0; i < len; i++) {
		double realPart = mersenne() / (double) mersenne.max();
		double imagPart = mersenne() / (double) mersenne.max();
		states[i] = complexd(realPart, imagPart);
		norm += pow(realPart, 2) + pow(imagPart, 2);
	}
	// нормирование вектора
	double norms[size];
	MPI_Gather(&norm, 1, MPI_DOUBLE, norms, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (!rank) {
		norm = 0;
		for (int i = 0; i < size; i++)
			norm += norms[i];
		norm = sqrt(norm);
	}
	MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for (long long i = 0; i < len; i++)
		states[i] /= norm;
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

// проверяет на 0 или 1 бит на numQubit-ном месте числа number
int checkBit(unsigned long long number, int countQubits, int numQubit) {
	unsigned long long mask = 1 << (countQubits - numQubit);
	if ((number & mask) == 0)
		return 0;
	return 1;
}

// возвращает число в 10-ой с.с. , которое получится после постановки 0 или 1 на numQubit-ный бит числа number
unsigned long long putZeroOrOne(int zeroOrOne, unsigned long long number, int countQubits, int numQubit) {
	unsigned long long mask = 1 << (countQubits - numQubit);
	if (zeroOrOne == 1)
		return (number | mask);
	mask = ~mask;
	return (number & mask);
}

// функция однокубитного преобразования
complexd *oneQubitTransformation(complexd *states, complexd *statesFriend, int countQubits, int numQubit, int rank, int rankFriend) {
	complexd *transformedStates = new complexd[len];
	vector <complexd> matrix[2] {{M_SQRT1_2, M_SQRT1_2}, {M_SQRT1_2, - M_SQRT1_2}};
	int iDiff = len * rank;
	int iDiffFriend = len * rankFriend;
	for (long long i = 0; i < len; i++) {
		unsigned long long index0 = putZeroOrOne(0, iDiff + i, countQubits, numQubit);
		complexd elem0 = (index0 >= iDiff && index0 < len + iDiff) ? states[index0 - iDiff] : statesFriend[index0 - iDiffFriend];
		unsigned long long index1 = putZeroOrOne(1, iDiff + i, countQubits, numQubit);
		complexd elem1 = (index1 >= iDiff && index1 < len + iDiff) ? states[index1 - iDiff] : statesFriend[index1 - iDiffFriend];
		transformedStates[i] = matrix[checkBit(iDiff + i, countQubits, numQubit)][0] * elem0 +
			matrix[checkBit(iDiff + i, countQubits, numQubit)][1] * elem1;
	}
	return transformedStates;
}

void testing(complexd *firstStates, complexd *secondStates, int size, int rank) {
	// сравнение векторов
	bool result = true;
	for (int i = 0; i < len; i++)
		if (firstStates[i] != secondStates[i]) {
			result = false;
			break;
		}
	// отсылка результатов сравнения со всех процессов 0-ому
	bool results[size];
	MPI_Gather(&result, 1, MPI_BYTE, results, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
	if (!rank) {
		for (int i = 0; i < size; i++)
			if (results[i] == false) {
				result = false;
				break;
			}
		if (result == true)
			cout << "✅ " << size << " threads" << endl;
		else
			cout << "❌ " << size << " threads" << endl;
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

	int countQubits;

	if (!strcmp(argv[1], "read") || !strcmp(argv[1], "test")) { // режимы чтения и тестирования
		// обработка ошибок с параметрами
		switch (argc) {
			case 2: return handleError(3, rank, 0, 0);
			case 3: return handleError(4, rank, 0, 0);
			case 4: return handleError(6, rank, 0, 0);
			case 5:
				if (!strcmp(argv[1], "read"))
					return handleError(7, rank, 0, 0);
				else
					return handleError(9, rank, 0, 0);
		}
		// общее для двух режимов чтение вектора
		countQubits = atoi(argv[3]);
		len = pow(2, countQubits) / size; // потому что кол-во процессов - степень 2
		complexd *states = new complexd[len];
		int numQubit = atoi(argv[4]);
		// generateVector(states, size, rank); // генерация вектора для тестирования времени
		if (readVectorFromFile(states, argv[2], size, rank))
			return handleError(5, rank, &states, 0);
		// подсчет потока для обмена данными
		int rankFriend;
		if (numQubit > log2(size))
			rankFriend = rank;
		else if (rank % (int)(size / pow(2, numQubit - 1)) < size / pow(2, numQubit))
			rankFriend = rank + size / pow(2, numQubit);
		else
			rankFriend = rank - size / pow(2, numQubit);
		// обмен данными с нужными потоком
		complexd *statesFriend = new complexd[len];
		for (int i = 0; i < len; i++)
			statesFriend[i] = states[i];
		if (rankFriend != rank) {
			MPI_Isend(states, len, MPI_DOUBLE_COMPLEX, rankFriend, 0, MPI_COMM_WORLD, &request);
			MPI_Recv(statesFriend, len, MPI_DOUBLE_COMPLEX, rankFriend, 0, MPI_COMM_WORLD, &status);
		}
		// общее для двух режимов преобразование вектора
		double timeStart = MPI_Wtime();
		complexd *transformedStates = oneQubitTransformation(states, statesFriend, countQubits, numQubit, rank, rankFriend);
		double timeEnd = MPI_Wtime();
		double time = timeEnd - timeStart;
		double maxTime;
		// вывод времени работы в файл для тестирования времени
		// MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		// if (!rank) {
		// 	ofstream fout(argv[5], ios_base::app);
		// 	fout << "\t" << numQubit << "\t" << countQubits << "\t" << size << "\t" << maxTime << endl;
		// 	fout.close();
		// }
		if (!strcmp(argv[1], "read")) { // вывод преобразованного вектора для режима чтения
			if (writeVectorInFile(transformedStates, argv[5], size, rank)) {
				delete [] transformedStates;
				return handleError(8, rank, &states, &statesFriend);
			}
		} else { // сравнение преобразованного вектора с вектором из файла в режиме тестирования
			complexd *validStates = new complexd[len];
			if (readVectorFromFile(validStates, argv[5], size, rank)) {
				delete [] transformedStates;
				delete [] validStates;
				return handleError(10, rank, &states, &statesFriend);
			}
			testing(transformedStates, validStates, size, rank);
			delete [] validStates;
		}
		delete [] transformedStates;
		delete [] states;
		delete [] statesFriend;
	} else if (!strcmp(argv[1], "generate")) { // режим генерации вектора в файл
		switch (argc) {
			case 2: return handleError(7, rank, 0, 0);
			case 3: return handleError(4, rank, 0, 0);
		}
		countQubits = atoi(argv[3]);
		len = pow(2, countQubits) / size;
		complexd *states = new complexd[len];
		generateVector(states, size, rank);
		if (writeVectorInFile(states, argv[2], size, rank))
			return handleError(8, rank, &states, 0);
		delete [] states;
	} else {
		return handleError(2, rank, 0, 0);
	}
	MPI_Finalize();
	return 0;
}
