#include <iostream>
#include "Matrix.h"

int main()
{
	Matrix matr( "inputData.csv" );
	matr.prepare();
    matr.reduceAsColumns();

	return EXIT_SUCCESS;
}
