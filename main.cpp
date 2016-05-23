#include <iostream>
#include "Matrix.h"

int main()
{
	Matrix matr( "inputData.csv" );
    matr.reduceAll();

    std::cerr << "\nMatrix after matr.reduceAll(); :\n";

    matr.printMatrix( matr.getPreparedMatrix(), std::cerr );


	return EXIT_SUCCESS;
}
