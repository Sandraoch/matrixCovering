#include <iostream>
#include "Matrix.h"

int main()
{
	Matrix matr( "inputData.csv" );
	matr.prepare();
	return EXIT_SUCCESS;
}
