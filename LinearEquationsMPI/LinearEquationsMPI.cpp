// LinearEquationsMPI.cpp : Решение СЛАУ методом Якоби (простых итераций) на MPI
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include "MatrixGeneration.h"
#include "mpi.h"
using namespace std;
int msgtag = 26;
MPI_Status status;
#define comm MPI_COMM_WORLD
int rankp, sizep;
vector<int> sizeOfMatrices = { 100, 200, 500 };
const int N = sizeOfMatrices[0];
double eps = 1e-5;
template<typename T>
void PrintArray(T** arr, int countOfRows, string fileName);
template<typename T>
void PrintArray(T* arr, int countOfRows, string fileName);
template<typename T>
void GetRowsFromMatrix(T** newArr, const vector<vector<T>>& commonArr, int countOfRows, int distance);
int main(int argc, char** argv)
{
	#pragma region HEADER
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        return 1;
    }
    if (MPI_Comm_size(MPI_COMM_WORLD, &sizep) != MPI_SUCCESS)
    {
        MPI_Finalize();
        return 2;
    }
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rankp) != MPI_SUCCESS)
    {
        MPI_Finalize();
        return 3;
    }
#pragma endregion


	#pragma region Считывание матриц из файла в векторы и конвертация в массивы
	//if (!MatrixGeneration()) // Matrix Generation
	//{
	//	cout << "Files were created" << endl;
	//}
		string
		fileNameA = "matrixA" + to_string(N) + ".txt", //a file name of matrices
		fileNameB = "matrixB" + to_string(N) + ".txt";
	double sum = 0.0, xMax = -1;
	ifstream matrixA, matrixB;
	matrixA.open(fileNameA); //opening of files
	matrixB.open(fileNameB);
	string str;  //for parsing data
	vector<vector<double>> a(N); // a[n][n]
	vector<double> b(N); // b[n]
	bool correctOpening = true;
	// Считывание из файлов матрицы А и вектора B
	if (matrixA.is_open())
	{
		int i = 0;
		while (getline(matrixA, str))
		{
			stringstream stream(str);
			while (getline(stream, str, '\t'))
			{
				a[i].push_back(stod(str));
			}
			i++;

		};
		matrixA.close();
	}
	else
	{
		cout << "A Matrix A is absent " << endl;
		correctOpening = false;
	}
	if (matrixB.is_open())
	{
		int i = 0;
		while (getline(matrixB, str))
		{
			stringstream stream(str);
			while (getline(stream, str, '\t'))
			{
				b[i++] = stod(str);
			}
		};
		matrixB.close();
	}
	else
	{
		cout << "A Matrix B is absent" << endl;
		correctOpening = false;
	}
#pragma endregion

	#pragma region Получение строк из общей матрицы для каждого процесса
	//Поделим строки матрицы на каждый процесс
	int l1 = N / sizep, l2 = N % sizep;
	int* kol = new int[sizep];
	for (int i = 0; i < sizep; i++)
	{
		kol[i] = l1;
	}
	if (l2) {
		if (rankp < l2) l1++;
		for (int i = 0; i < l2; i++) kol[i]++;
	}
	int countOfRows = l1;
	int* displacement = new int[sizep] {0};
	for (int i = 1; i < sizep; i++)
	{
		displacement[i] = displacement[i - 1] + kol[i - 1];
	}
	int distance = 0; // переменная, чтобы узнать, с какой строчки начинать брать каждому процессу от общей матрицы
	for (int i = 0; i < rankp; i++)
	{
		distance += kol[i]; // каждый процесс возьмёт только свои строчки
	}
	double** processArray = new double* [N];
	for (int i = 0; i < countOfRows; i++)
	{
		processArray[i] = new double[N];
	}
	GetRowsFromMatrix(processArray, a, countOfRows, distance); 
#pragma endregion	
															
	#pragma region Расчёты
	double* xNew = new double[N] {0};
	double* xOld = new double[N] {0};
	double globalDeviation = 0.0;
	if (correctOpening)
	{		
		do
		{
			xMax = -1;
			for (int i = 0; i < countOfRows; i++) {
				sum = 0;
				for (int j = 0; j < N; j++)
				{
					if (i + distance != j)
					{
						sum += processArray[i][j] * xNew[j];
					}
				}
				xOld[i + distance] = 1.0 / processArray[i][i + distance] * (b[i + distance] - sum);

				//нахождение максимального отклонения
				if (fabs(xOld[i + distance] - xNew[i + distance]) > xMax)
				{
					xMax = fabs(xOld[i + distance] - xNew[i + distance]);
				}
			}
			//синхронизация всех процессов
			MPI_Barrier(comm);
			// Нахождение глобального отклонения
			MPI_Allreduce(&xMax, &globalDeviation, 1, MPI_DOUBLE, MPI_MAX, comm);
			MPI_Allgatherv(&xOld[distance], kol[rankp], MPI_DOUBLE, xNew, kol, displacement, MPI_DOUBLE, comm);
		} while (globalDeviation > eps);
	}
	else
	{
		cout << "Files wasn't found" << endl;
	}
	#pragma endregion

	if (!rankp)
	{		
		PrintArray(xNew, N, "result" + to_string(N) + ".txt");
		cout << "Program's worked" << endl;
	}
	MPI_Barrier(comm);
	#pragma region Footer
	for (int i = 0; i < countOfRows; i++)
	{
		delete[] processArray[i];
	}
	delete[] xNew, xOld, displacement, kol;
	MPI_Finalize();	
	
	return MPI_SUCCESS;
}
template<typename T>
void PrintArray(T** arr, int countOfRows, string fileName)
{
	ofstream output(fileName);
	if (output.is_open())
	{
		for (int i = 0; i < countOfRows; i++)
		{
			bool sep = false;
			for (int j = 0; j < N; j++)
			{
				if (sep)
				{
					output << " ";
				}
				sep = true;
				output << setw(10) << arr[i][j];
			}
			output << endl;
		}
		output.close();
	}
	else
	{
		cout << "File is absent" << endl;
	}

}
template<typename T>
void PrintArray(T* arr, int countOfRows, string fileName)
{
	ofstream output(fileName);
	if (output.is_open())
	{
		for (int i = 0; i < countOfRows; i++)
		{
			output << i + 1 << ". " << setw(15) << arr[i] << endl;
		}
		output.close();
	}
	else
	{
		cout << "File is absent" << endl;
	}

}
template<typename T>
void GetRowsFromMatrix(T** newArr, const vector<vector<T>>& commonArr, int countOfRows, int distance)
{
	for (int i = 0; i < countOfRows; i++)
	{
		for (int j = 0; j < N; j++)
		{
			newArr[i][j] = commonArr[distance][j];
		}
		distance++;
	}
}