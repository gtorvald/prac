#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <string>
#include <omp.h>

using namespace std;

// количество кубитов
#define N 4
// номер кубита , по которому проводится преобразование
#define NUM_QUBIT 1

// возвращает двоичное представление числа num длиной len
string getBinary(int number, int len) {
	string result;
	while (number > 0) {
		if (number % 2 == 0)
			result = '0' + result;
		else
			result = '1' + result;
		number /= 2;
		len--;
	}
	while (len--)
		result = '0' + result;
	return result;
}

// возвращает число в 10-ой с.с. , которое получится после постановки 0 или 1 
// на numQubit-ное место в 2-ое представление числа number
int putZeroOrOne(char zeroOrOne, string &number, int numQubit) {
	number[numQubit - 1] = zeroOrOne;
	int numberNew = 0;
	for (int i = 0; i < number.length(); i++)
		if (number[i] == '1')
			numberNew += pow(2, number.length() - i - 1);
	return numberNew;
}

// функция однокубитного преобразования
vector <complex <double> > oneQubitTransformation(vector <complex <double> > &vec, vector <vector <complex <double> > > &matrix, int numQubit) {
	vector <complex <double> > vecTransformed(vec.size());
	#pragma omp parallel for
	for (int i = 0; i < vec.size(); i++) {
		string iTwoBit = getBinary(i, log(vec.size()) / log(2));
		vecTransformed[i] = matrix[iTwoBit[numQubit] - '0'][0] * vec[putZeroOrOne('0', iTwoBit, numQubit)] +
			matrix[iTwoBit[numQubit] - '0'][1] * vec[putZeroOrOne('1', iTwoBit, numQubit)]; 
	}
	return vecTransformed;
}

int	main() {
	// инициализация вектора случайными значениями
	int vectorSize = pow(2, N);
	vector <complex <double> > vec(vectorSize);
	#pragma omp parallel for
	for (int i = 0; i < vec.size(); i++)
		vec[i] = complex<double> (1. / rand(), 1. / rand());
	// инициализация матрицы преобразования Адамара
	vector <vector <complex <double> > > matrix(2);
	matrix[0].resize(2);
	matrix[1].resize(2);
	#pragma omp parallel for collapse(2)
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++) {
			matrix[i][j] = 1. / sqrt(2);
			if (i == 1 && j == 1)
				matrix[i][j] = - matrix[i][j];
		}
	// вызов функции однокубитного преобразования
	vector <complex <double> > vecTransformed = oneQubitTransformation(vec, matrix, NUM_QUBIT);
	// вывод нового вектора
	for (int i = 0; i < vecTransformed.size(); i++)
		cout << vecTransformed[i] << endl;
	return 0;
}