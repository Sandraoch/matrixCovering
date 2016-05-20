#pragma once
#ifndef MATRIX_H
#define MATRIX_H
#include <vector>

class Ellement_t {
public:
	size_t row, col;
	unsigned char value;
	Ellement_t(size_t row, size_t col, unsigned char val);
};

class Matrix {
public:
	typedef std::vector<std::vector<Ellement_t>> Matrix_t;
private:
	Matrix_t mOriginal;
	Matrix_t mPrepared;
	std::vector<size_t> preparedRows;
    void deleteRowsAndCoveredColumns(std::vector<size_t>& r, Matrix_t &m);
    std::vector<size_t> getKernelRows( const Matrix_t& m );
    std::vector<size_t> getReducingColumns( const Matrix_t& m );
    void deleteColumns(std::vector<size_t>& c, Matrix_t &m);
public:
    Matrix_t getOriginalMatrix();
    Matrix_t getPreparedMatrix();
	void printMatrix( const Matrix_t & m, std::ostream & oStr );
	void prepare();
    void reduceAsColumns();
	Matrix(std::string filePath);
};

#endif