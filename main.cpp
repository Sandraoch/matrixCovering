#include <iostream>
//#define DEBUG_MODE
#include "Matrix.h"
#include <fstream>

int main()
{
    //Matrix matr("inputData_with_problem.csv");
    Matrix matr("inputData.csv");

    matr.reduceAll();

    std::fstream fStream( "res.txt" );

    Matrix::Matrix_t m = matr.getOriginalMatrix();

    fStream << "\nOriginal matrix:\n";
    matr.printMatrix(m, fStream);

    fStream << "\nCovering is:\n";

    Matrix::ListOfIndexes_t covering = matr.getCurCovering();

    for( const auto &row : covering )
    {
        fStream << int(m[row][0].row + 1) << "|\t";
        for( const auto &el : m[row] )
            fStream << (int)el.value  << '\t';
        fStream << std::endl;
    }

    fStream << "\nGradient method without preparing:\n";
    matr.makePreparedMatrixOriginal();
    matr.gradientMethod();
    covering = matr.getCurCovering();
    for( const auto &row : covering )
    {
        fStream << int(m[row][0].row + 1) << "|\t";
        for( const auto &el : m[row] )
            fStream << (int)el.value  << '\t';
        fStream << std::endl;
    }

    fStream << "\nGradient method with preparing:\n";
    matr.makePreparedMatrixOriginal();
    matr.prepare();
    matr.reduceAsColumns();
    matr.gradientMethod();
    covering = matr.getCurCovering();
    for( const auto &row : covering )
    {
        fStream << int(m[row][0].row + 1) << "|\t";
        for( const auto &el : m[row] )
            fStream << (int)el.value  << '\t';
        fStream << std::endl;
    }

    return EXIT_SUCCESS;
}
