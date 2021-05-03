#include <iostream>
#include <vector>
#include <complex>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include <omp.h>

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
			case 6: cout << "Error 06: invalid parameters: no eps" << endl; break;
			case 7: cout << "Error 07: invalid parameters: no outputFile" << endl; break;
			case 8: cout << "Error 08: invalid outputFile" << endl; break;
			case 9: cout << "Error 09: invalid parameters: no file with valid states" << endl; break;
			case 10: cout << "Error 10: invalid file with valid states" << endl; break;
			case 11: cout << "Error 11: invalid parameters: no timeFile" << endl; break;
			case 12: cout << "Error 12: invalid parameters: no count of threads for OpenMP" << endl; break;
		}
		cout << "usage:" << endl;
		cout << "\tmpirun -np <threads> go read <inputFile> <countQubits> <eps> <outputFile> <timeFile> <countThreads>" << endl;
		cout << "\tmpirun -np <threads> go generate <outputFile> <countQubits> <countThreads>" << endl;
		cout << "\tmpirun -np <threads> go test <inputFile> <countQubits> <eps> <validFile> <countThreads>" << endl;
	}
	MPI_Finalize();
	if (states)
		delete [] *states;
	if (statesFriend)
		delete [] *statesFriend;
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

// функция однокубитного преобразования
complexd *oneQubitTransformation(complexd *states, complexd matrix[2][2], complexd *statesFriend, int countQubits, int numQubit, int rank, int rankFriend) {
	complexd *transformedStates = new complexd[len];
	int diff = countQubits - numQubit;
	unsigned long long mask = 1 << diff;

	if (rank == rankFriend)
		for (long long i = 0; i < len; i++) {
			int index = (i & mask) >> diff;
			transformedStates[i] = matrix[index][0] * states[i & ~mask] + matrix[index][1] * states[i | mask];
		}
	if (rank < rankFriend)
		for (long long i = 0; i < len; i++)
			transformedStates[i] = matrix[0][0] * states[i] + matrix[0][1] * statesFriend[i];
	else if (rank > rankFriend)
		for (long long i = 0; i < len; i++)
			transformedStates[i] = matrix[1][0] * statesFriend[i] + matrix[1][1] * states[i];

	return transformedStates;
}

double testing(complexd *firstStates, complexd *secondStates, int size, int rank) {
	// сравнение векторов
	complexd localSum = complexd(0, 0);
	for (int i = 0; i < len; i++)
		localSum += firstStates[i] * conj(secondStates[i]);
	// отсылка результатов сравнения со всех процессов 0-ому
	complexd sum = 0;
	MPI_Reduce(&localSum, &sum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
	return 1 - pow(pow(sum.real(), 2) + pow(sum.imag(), 2), 2);
}

double normal_dis_gen() {
	double S = 0;
	for (int i = 0; i < 12; ++i) {
		S += (double) rand() / RAND_MAX;
	}
	return S - 6;
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
			case 5: return (!strcmp(argv[1], "read")) ? handleError(7, rank, 0, 0) : handleError(9, rank, 0, 0);
			case 6: return (!strcmp(argv[1], "read")) ? handleError(11, rank, 0, 0) : handleError(12, rank, 0, 0);
			case 7:
				if (!strcmp(argv[1], "read"))
					return handleError(12, rank, 0, 0);
		}
		omp_set_num_threads(atoi(argv[7]));
		// общее для двух режимов чтение вектора
		countQubits = atoi(argv[3]);
		len = pow((double) 2, countQubits) / size; // потому что кол-во процессов - степень 2
		complexd *states = new complexd[len];
		// generateVector(states, size, rank); // генерация вектора для тестирования времени
		if (readVectorFromFile(states, argv[2], size, rank))
			return handleError(5, rank, &states, 0);
		// преобразование n-Адамар
		complexd *transformedStates;
		complexd *statesFriend;
		complexd matrix[2][2];
		matrix[0][0] = M_SQRT1_2;
		matrix[0][1] = M_SQRT1_2;
		matrix[1][0] = M_SQRT1_2;
		matrix[1][1] = - M_SQRT1_2;
		// в режите чтения преобразование запускается с шумом
		// в тестовом для расчета потери точности - без шума
		if (!strcmp(argv[1], "read")) {
			double eps = atof(argv[4]);
			double phi;
			if (rank == 0)
				phi = normal_dis_gen();
			MPI_Bcast(&phi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			matrix[0][0] = M_SQRT1_2 * (cos(eps * phi) - sin(eps * phi));
			matrix[0][1] = M_SQRT1_2 * (cos(eps * phi) + sin(eps * phi));
			matrix[1][0] = matrix[0][1];
			matrix[1][1] = M_SQRT1_2 * (sin(eps * phi) - cos(eps * phi));
		}
		double timeStart = MPI_Wtime();
		for (int numQubit = 1; numQubit <= countQubits; numQubit++) {
			// подсчет потока для обмена данными
			int rankFriend;
			long double sizePow2numQubit1 = size / pow((double) 2, (int) numQubit - 1);
			long double sizePow2numQubit = sizePow2numQubit1 / 2;
			if (numQubit > log2(size))
				rankFriend = rank;
			else if (rank % (int) sizePow2numQubit1 < sizePow2numQubit)
				rankFriend = rank + sizePow2numQubit;
			else
				rankFriend = rank - sizePow2numQubit;
			if (rankFriend != rank) {
				// обмен данными с нужными потоком
				MPI_Isend(states, len, MPI_DOUBLE_COMPLEX, rankFriend, 0, MPI_COMM_WORLD, &request);
				MPI_Recv(statesFriend, len, MPI_DOUBLE_COMPLEX, rankFriend, 0, MPI_COMM_WORLD, &status);
				// общее для двух режимов преобразование вектора
				transformedStates = oneQubitTransformation(states, matrix, statesFriend, countQubits, numQubit, rank, rankFriend);
			}
			else
				// общее для двух режимов преобразование вектора
				transformedStates = oneQubitTransformation(states, matrix, states, countQubits, numQubit, rank, rankFriend);
			delete [] states;
			states = transformedStates;
			transformedStates = NULL;
		}
		double timeEnd = MPI_Wtime();
		double time = timeEnd - timeStart;
		double maxTime;
		delete [] statesFriend;
		// вывод времени работы в файл для тестирования времени
		// MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		// if (rank == 0) {
		// 	ofstream fout(argv[6], ios_base::app);
		// 	fout << size << ' ' << maxTime << endl;
		// 	fout.close();
		// }
		if (!strcmp(argv[1], "read")) { // вывод преобразованного вектора для режима чтения
			if (writeVectorInFile(states, argv[5], size, rank))
				return handleError(8, rank, &states, 0);
		} else { // сравнение преобразованного вектора с вектором из файла в режиме тестирования
			complexd *validStates = new complexd[len];
			if (readVectorFromFile(validStates, argv[5], size, rank)) {
				delete [] validStates;
				return handleError(10, rank, &states, 0);
			}
			double precision = testing(states, validStates, size, rank);
			if (rank == 0) {
				ofstream fout(argv[7], ios_base::app);
				fout << precision << endl;
				fout.close();
			}
			delete [] validStates;
		}
		delete [] states;
		delete [] statesWithoutNoise;
	 } else if (!strcmp(argv[1], "generate")) { // режим генерации вектора в файл
		switch (argc) {
			case 2: return handleError(7, rank, 0, 0);
			case 3: return handleError(4, rank, 0, 0);
		}
		countQubits = atoi(argv[3]);
		len = pow((double) 2, countQubits) / size;
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
