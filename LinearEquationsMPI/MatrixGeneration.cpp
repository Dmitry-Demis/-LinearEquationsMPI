#include "MatrixGeneration.h"

bool MatrixGeneration()
{
	srand((unsigned)time(NULL));
	string fileNameA, fileNameB;
	vector<int> matrixSize = {100, 200, 500};
	for (const auto& N : matrixSize)
	{
		fileNameA = "matrixA" + to_string(N) + ".txt";
		fileNameB = "matrixB" + to_string(N) + ".txt";
		ofstream matrixA(fileNameA), matrixB(fileNameB);
		if (!matrixA.is_open() || !matrixB.is_open())
		{
			cout << "File recording error" << endl;
			return false;
		}
		vector<double> matrA(N), matrB(N);
		bool sep = false;
		for (int i = 0; i < N; i++)
		{
			double sum = 0.0;
			for (int j = 0; j < N; j++)
			{
				matrA[j] = rand() % N + 0.5;
				sum += matrA[j];
			}
			matrA[i] = sum + 1;
			for (int j = 0; j < N; j++)
			{
				matrixA << setw(4) << matrA[j] << "\t";
			}
			matrixA << endl;
			matrB[i] = rand() % N + 0.7;

			if (sep)
			{
				matrixB << "\t";
			}
			sep = true;
			matrixB << matrB[i];
		}
		matrixA.close();
		matrixB.close();
		matrA.clear();
		matrB.clear();
	}
	return true;
}