#include <iostream>
#include <complex>
#include <mpi.h>

using namespace std;

typedef complex<double> complexd;

complexd *oneQubitTransformation(complexd *states, unsigned long long len, complexd matrix[2][2], int countQubits, int numQubit, int rank, int size) {
	MPI_Status status;
	MPI_Request request;
	complexd *statesFriend = new complexd[len];
	int rankFriend;
	complexd *transformedStates = new complexd[len];
	int diff = countQubits - numQubit;
	unsigned long long mask = 1L << diff;
	if ((((unsigned long long) rank * len) & mask) == 0)
		rankFriend = (((unsigned long long) rank * len) | mask) / len;
	else
		rankFriend = (((unsigned long long) rank * len) xor mask) / len;
	if (rankFriend != rank) {
		MPI_Isend(states, len, MPI_DOUBLE_COMPLEX, rankFriend, 0, MPI_COMM_WORLD, &request);
		MPI_Recv(statesFriend, len, MPI_DOUBLE_COMPLEX, rankFriend, 0, MPI_COMM_WORLD, &status);
	}

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

	delete [] statesFriend;
	return transformedStates;
}

void H(complexd **states, int countQubits, int numQubit, int rank, int size) {
	unsigned long long len = pow((double) 2, countQubits) / size;
	complexd matrix[2][2];
	matrix[0][0] = M_SQRT1_2;
	matrix[0][1] = M_SQRT1_2;
	matrix[1][0] = M_SQRT1_2;
	matrix[1][1] = - M_SQRT1_2;

	complexd *transformedStates = oneQubitTransformation(*states, len, matrix, countQubits, numQubit, rank, size);
	delete [] *states;
	*states = transformedStates;
}

void Hn(complexd **states, int countQubits, int rank, int size) {
	for (int i = 1; i <= countQubits; i++)
		H(states, countQubits, i, rank, size);
}

void ROT(complexd **states, int countQubits, int numQubit, int rank, int size) {
	unsigned long long len = pow((double) 2, countQubits) / size;
	unsigned long long mask = 1LL << (countQubits - numQubit); 
	unsigned long long diff = rank * len;

	for (int i = 0; i < len; i++)
		if (((i + diff) & mask) != 0)
			(*states)[i] = - (*states)[i];
}

void CROT(complexd **states, int countQubits, int firstQubit, int secondQubit, int rank, int size) {
	unsigned long long len = pow((double) 2, countQubits) / size;
	unsigned long long mask = (1LL << (countQubits - firstQubit)) | (1LL << (countQubits - secondQubit)); 
	unsigned long long diff = rank * len;

	for (int i = 0; i < len; i++)
		if (((i + diff) & mask) == mask)
			(*states)[i] = - (*states)[i];
}

void NOT(complexd **states, int countQubits, int numQubit, int rank, int size) {
	unsigned long long len = pow((double) 2, countQubits) / size;
	complexd matrix[2][2];
	matrix[0][0] = 0;
	matrix[0][1] = 1;
	matrix[1][0] = 1;
	matrix[1][1] = 0;

	complexd *transformedStates = oneQubitTransformation(*states, len, matrix, countQubits, numQubit, rank, size);
	delete [] *states;
	*states = transformedStates;
}

void CNOT(complexd **states, int countQubits, int firstQubit, int secondQubit, int rank, int size) {
	MPI_Status status;
	MPI_Request request;
	unsigned long long len = pow((double) 2, countQubits) / size;
	unsigned long long mask = (1LL << (countQubits - firstQubit)) | (1LL << (countQubits - secondQubit)); 
	unsigned long long diff = rank * len;
	unsigned long long rankFriend;
	complexd *statesFriend = new complexd[len];
	complexd *transformedStates = new complexd[len];

	if ((diff & mask) == (1LL << (countQubits - firstQubit)))
		rankFriend = (diff | mask) / len;
	else if ((diff & mask) == mask)
		rankFriend = (diff  xor (1LL << (countQubits - secondQubit))) / len;
	else
		rankFriend = rank;
	if (rankFriend != rank) {
		MPI_Isend(*states, len, MPI_DOUBLE_COMPLEX, rankFriend, 0, MPI_COMM_WORLD, &request);
		MPI_Recv(statesFriend, len, MPI_DOUBLE_COMPLEX, rankFriend, 0, MPI_COMM_WORLD, &status);
	}
	unsigned long long diffFriend = rankFriend * len;

	for (int i = 0; i < len; i++)
		if (((i + diff) & mask) == (1LL << (countQubits - firstQubit))) {
			if (rank != rankFriend)
				transformedStates[i] = statesFriend[i];
			else
				transformedStates[i] = (*states)[((i + diff) | mask) - diff];
		}
		else if (((i + diff) & mask) == mask) {
			if (rank != rankFriend)
				transformedStates[i] = statesFriend[i];
			else
				transformedStates[i] = (*states)[((i + diff) xor (1LL << (countQubits - secondQubit))) - diff];
		}
		else
			transformedStates[i] = (*states)[i];
	delete [] statesFriend;
	delete [] *states;
	*states = transformedStates;
}
