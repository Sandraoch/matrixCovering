#include <iostream>
//#define DEBUG_MODE
#include "Matrix.h"
#include <fstream>

int main()
{
    //Matrix matr("inputData_with_problem.csv");
    Matrix matr("inputData.csv");

    matr.reduceAll();

    Matrix::Matrix_t m = matr.getOriginalMatrix();

    matr.printMatrix(m, std::cerr);

    std::cerr << "\nCovering is:\n";

    Matrix::ListOfIndexes_t covering = matr.getCurCovering();

    for( const auto &row : covering )
    {
        std::cerr << int(m[row][0].row + 1) << "|\t";
        for( const auto &el : m[row] )
            std::cerr << (int)el.value  << '\t';
        std::cerr << std::endl;
    }

    return EXIT_SUCCESS;
}
