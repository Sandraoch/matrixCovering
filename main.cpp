#include <iostream>
//#define DEBUG_MODE
#include "Matrix.h"
#include <fstream>

int main()
{
	Matrix matr("inputData.csv");
	matr.reduceAll();

	std::cerr << "\nMatrix after matr.reduceAll(); :\n";

	std::ofstream out("res.txt");

	matr.printMatrix(matr.getPreparedMatrix(), out);

	return EXIT_SUCCESS;
}
