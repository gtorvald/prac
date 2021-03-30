#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <string>
#include <omp.h>
#include <fstream>

using namespace std;

// количество кубитов
#define COUNT_QUBIT 20
// номер кубита , по которому проводится преобразование
#define NUM_QUBIT 1

// проверяет на 0 или 1 бит на numQubit-ном месте числа number
int checkBit(unsigned long long number, int numQubit) {
	unsigned long long mask = 1 << (COUNT_QUBIT - numQubit);
	if ((number & mask) == 0)
		return 0;
	return 1;
}

// возвращает число в 10-ой с.с. , которое получится после постановки 0 или 1 на numQubit-ный бит числа number
unsigned long long putZeroOrOne(int zeroOrOne, unsigned long long number, int numQubit) {
	unsigned long long mask = 1 << (COUNT_QUBIT - numQubit);
	if (zeroOrOne == 1)
		return (number | mask);
	mask = ~mask;
	return (number & mask);
}

// функция однокубитного преобразования
vector <complex <double> > oneQubitTransformation(vector <complex <double> > &vec, vector <vector <complex <double> > > &matrix, int numQubit) {
	vector <complex <double> > vecTransformed(vec.size());
	#pragma omp parallel for
	for (long long i = 0; i < vec.size(); i++)
		vecTransformed[i] = matrix[checkBit(i, numQubit)][0] * vec[putZeroOrOne(0, i, numQubit)] +
			matrix[checkBit(i, numQubit)][1] * vec[putZeroOrOne(1, i, numQubit)]; 
	return vecTransformed;
}

int	main(int argc, char **argv) {
	// инициализация вектора случайными значениями c последующим нормированием
	unsigned long long vectorSize = pow(2, COUNT_QUBIT);
	vector <complex <double> > vec(vectorSize);
	double norm = 0;
	#pragma omp parallel
	{
		#pragma omp for reduction(+:norm)
		for (long long i = 0; i < vec.size(); i++) {
			double realPart = 1. / rand();
			double imagPart = 1. / rand();
			vec[i] = complex<double> (realPart, imagPart);
			norm += pow(realPart, 2) + pow(imagPart, 2);
		}
		norm = sqrt(norm);
		#pragma omp for
		for (long long i = 0; i < vec.size(); i++)
			vec[i] = vec[i] / norm;
	}
	// инициализация матрицы преобразования Адамара
	vector <vector <complex <double> > > matrix(2);
	matrix[0].resize(2);
	matrix[1].resize(2);
	#pragma omp parallel for collapse(2)
	for (long long i = 0; i < 2; i++)
		for (long long j = 0; j < 2; j++) {
			matrix[i][j] = 1. / sqrt(2);
			if (i == 1 && j == 1)
				matrix[i][j] = - matrix[i][j];
		}
	// вызов функции однокубитного преобразования
	double begin = omp_get_wtime();
	vector <complex <double> > vecTransformed = oneQubitTransformation(vec, matrix, NUM_QUBIT);
	double end = omp_get_wtime();
	// вывод времени работы программы в файл
	ofstream fout(argv[1], ios_base::app);
	fout << end - begin << endl;
	fout.close();
	// вывод нового вектора
	// for (long long i = 0; i < vecTransformed.size(); i++)
	// 	cout << vecTransformed[i] << endl;
	return 0;
}
