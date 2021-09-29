// LinearEquationsMPI.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include "mpi.h"
#include "MatrixGeneration.h"

using namespace std;
int msgtag = 26;
MPI_Status status;
#define comm MPI_COMM_WORLD
int rankp, sizep;
const int N = 11;
double eps = 1e-4;

template <typename T>
void ConvertVectorToArray(vector<vector<T>> vec, T** arr)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			arr[i][j] = vec[i][j];
		}
	}
}
template<typename T>
void ConvertVectorToArray(vector<T> vec, T* arr)
{
	for (int i = 0; i < N; i++)
	{
		arr[i] = vec[i];
	}
}
template<typename T>
void PrintArray(T** arr,  int countOfRows, string fileName)
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
			bool sep = false;

			if (sep)
			{
				output << " ";
			}
			sep = true;
			output << setw(10) << arr[i];
		}
		output.close();
	}
	else
	{
		cout << "File is absent" << endl;
	}

}

template<typename T>
void GetRowsFromMatrix(T** newArr, T** commonArr, int countOfRows, int distance)
{
	cout << "Process " << rankp << " has " << distance << "; " << distance + countOfRows << endl;
	for (int i = 0; i < countOfRows; i++)
	{
		for (int j = 0; j < N; j++)
		{
			newArr[i][j] = commonArr[distance][j];
		}
		distance++;
	}
}

#define P(x)  #x << " = " << x

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
	if (!MatrixGeneration()) // Matrix Generation
	{
		cout << "Files were created" << endl;
	}

	string
		fileNameA = "matrixA" + to_string(N) + ".txt", //a file name of matrices
		fileNameB = "matrixB" + to_string(N) + ".txt";
	double sum1 = 0.0, sum2 = 0.0, xMax = -1;
	ifstream matrixA, matrixB;
	matrixA.open(fileNameA); //opening of files
	matrixB.open(fileNameB);

	string str;  //for parsing data
	vector<vector<double>> a(N); // a[n][n]
	vector<double> b(N); // b[n]
	bool correctOpening = true, loading = true;
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

	//общая матрица А для всех процессов
	double** mA = new double* [N];
	for (int i = 0; i < N; i++)
		mA[i] = new double[N];
	double* mB = new double[N];
	ConvertVectorToArray(a, mA);
	ConvertVectorToArray(b, mB);
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
	int distance = 0;
	//int* displacements = new int[sizep] {0};
		for (int i = 0; i < rankp; i++)
		{
			//displacements[i+1] = distance;
			distance += kol[i]; // каждый процесс возьмёт только свои строчки
		}
	double** processArray = new double* [N];
	for (int i = 0; i < countOfRows; i++)
	{
		processArray[i] = new double[N];
	}
	GetRowsFromMatrix(processArray, mA, countOfRows, distance); // working 
#pragma endregion


	PrintArray(processArray, countOfRows, "output" + to_string(rankp) + fileNameA);
																
	#pragma region Расчёты

	#pragma endregion

	

	vector<double> matrixX(N, -1.0), matrixXX(N, -1.0);
	double* mX = new double[N];
	//	do{
			xMax = -1;
			int dist = distance;
			
			//if (rankp)
			{
				for (int i = 0; i < countOfRows; i++)
				{

					sum1 = 0;
					for (int j = 0; j < N; j++)
					{
						if (i + distance != j)
						{
							sum1 += processArray[i][j] * matrixX[j];
						}
					}					
					matrixXX[i + distance] = 1.0 / processArray[i][i + distance] * (mB[i + distance] - sum1);
					
					cout << P(i) << " " << P(rankp) << " " << P(xMax) << endl;
					if (fabs(matrixXX[i] - matrixX[i]) > xMax)
					{
						xMax = fabs(matrixXX[i] - matrixX[i]);
					}
					
				}				
				ConvertVectorToArray(matrixXX, mX);				
				MPI_Barrier(comm);
				// Нахождение глобального отклонения
				double globalDeviation = 0.0;			
				MPI_Allreduce(&xMax, &globalDeviation, 1, MPI_DOUBLE, MPI_MAX, comm);	
				// Всё выше работает



				double* result = new double[N] {0};
				const int countOfElements = countOfRows;

				// disp[size] = {...}
				// Allgatherv
				
			}

			matrixX = matrixXX;
		//} while (globalDeviation > eps);


    MPI_Finalize();
    return MPI_SUCCESS;
}

