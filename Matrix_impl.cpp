#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iostream>

#include "Matrix.h"

Ellement_t::Ellement_t(size_t row, size_t col, unsigned char val)
{
	this->row = row;
	this->col = col;
	this->value = val;
}

///Read matrix from .csv file with delims like "," and "\n"
Matrix::Matrix(std::string filePath)
{
	std::ifstream ifS( filePath );
	if( !ifS.is_open())
		throw std::runtime_error(filePath + " is not valid file ");
	
	std::string tempStrRow;

	size_t i = 0, j = 0;

	while (std::getline(ifS, tempStrRow, '\n'))
	{
		j = 0;
		std::stringstream sStr(tempStrRow);
		std::string tempStr;
		std::vector<Ellement_t> row;
		unsigned char tempU = 0;
		bool successfully = false;

		while (std::getline(sStr, tempStr, ','))
		{
            tempStr = tempStr[0];
			if (tempStr == "1")
				tempU = 1;
			else if (tempStr == "0")
				tempU = 0;
			else
                throw std::logic_error(
                        "Unknow symbol in file: " +
                        filePath +
                        ": '" +
                        tempStr + "'"
                        );

			Ellement_t newEll( i, j, tempU );

			row.push_back( newEll );
			++j;
			successfully = true;
		}

		if (successfully)
		{
			this->mOriginal.push_back(row);
			++i;
		}
	}
	this->mPrepared = this->mOriginal;
}

void Matrix::prepare()
{
	if (!mPrepared.size() || !mPrepared.at(0).size())
		throw std::logic_error("Empty matrix ");

    auto kernel = getKernelRows( mPrepared );
    preparedRows.insert( preparedRows.end(), kernel.begin(), kernel.end() );

#ifdef DEBUG_MODE
    std::cerr << "\nBefore deleteRowsAndCoveredColomns( preparedRows );\n";
    std::cerr << "preparedRows: " << std::endl;

    for( auto &col : preparedRows )
        std::cerr << (col + 1) << '\t';
    std::cerr << std::endl;
    printMatrix(mPrepared, std::cerr);
#endif //DEBUG_MODE

    deleteRowsAndCoveredColumns(preparedRows, mPrepared);

#ifdef DEBUG_MODE
    std::cerr << "\nAfter deleteRowsAndCoveredColomns(preparedRows);\n";
    printMatrix(mPrepared, std::cerr);
#endif //DEBUG_MODE

}

void Matrix::deleteRowsAndCoveredColumns( std::vector<size_t>& r, Matrix_t &m )
{
    size_t prepSize = m.size();
	std::sort(r.begin(), r.end());
	
	std::vector<size_t> deleteRows = r;
    std::vector<size_t> deleteColumns;
	
	for (const auto &row : r)
	{
        for (size_t j = 0; j < m[row].size(); ++j)
            if (m[row][j].value == 1)
                deleteColumns.push_back(j);
	}

	std::sort(deleteRows.begin(), deleteRows.end());
	for (int i = deleteRows.size() - 1; i >= 0; --i)
        std::swap(m[deleteRows[i]], m[prepSize - i - 1]);
	
    m.erase(m.begin() + (prepSize - deleteRows.size()), m.end());

    std::sort( deleteColumns.begin(), deleteColumns.end() );

#ifdef DEBUG_MODE
    std::cerr << "deleteColomns: " << std::endl;

    for( auto &col : deleteColumns )
        std::cerr << (col + 1) << '\t';
#endif
    this->deleteColumns( deleteColumns, mPrepared );
}

void Matrix::printMatrix(const Matrix_t & m, std::ostream & oStream )
{
	if (m.empty())
	{
		oStream << "\nEmpty matrix\n";
		return;
	}

	for (const auto & ell : m[0])
        oStream << '\t' << (ell.col + 1);

	for (const auto &row : m)
	{
		oStream << '\n';
        oStream << (row[0].row + 1) << "|\t";
		for (const auto & ell : row)
			oStream << (int)ell.value << '\t';
	}
    oStream << '\n';

}

void Matrix::deleteColumns( std::vector<size_t> &c, Matrix_t &m )
{
    size_t i = 0;
    for( const auto &col : c )
    {
        for( auto &row : m )
            row.erase( row.begin() + col - i );
        ++i;
    }
}

std::vector<size_t> Matrix::getKernelRows( const Matrix_t& m )
{
    std::vector<size_t> answ;
    for (size_t j = 0; j < m.at(0).size(); ++j)
    {
        size_t onesCount = 0;
        size_t curPreparedRow = -1;
        for (size_t i = 0; i < m.size(); i++)
        {
            if (m[i][j].value)
            {
                ++onesCount;
                curPreparedRow = i;
            }
        }

        if (onesCount == 1)
            answ.push_back(curPreparedRow);
    }

    return answ;
}

Matrix::Matrix_t Matrix::getOriginalMatrix()
{
    return mOriginal;
}

Matrix::Matrix_t Matrix::getPreparedMatrix()
{
    return mPrepared;
}

void Matrix::reduceAsColumns()
{
    std::vector<size_t> reducingColomns = getReducingColumns( mPrepared );

#ifdef DEBUG_MODE
    std::cerr << "\n Reducing colomns:\n";
    for( const auto &reduce : reducingColomns )
        std::cerr << (mPrepared[0][reduce].col + 1) << '\t';
    std::cerr << std::endl;
#endif //DEBUG_MODE
}

std::vector<size_t> Matrix::getReducingColumns( const Matrix_t &m )
{
    std::vector<size_t> answ;
    
    for( size_t i = 0; i < m.at(0).size(); ++i )
        for( size_t j = 0; j < m.at(0).size(); ++j )
        {
            bool covering = true;
            for( size_t k = 0; k < m.size(); ++k )
                if( m[k][j].value == 1 && m[k][i].value == 0 )
                {
                    covering = false;
                    break;
                }

            if( covering )
            {
                size_t onesOnFst = 0, onesInSnd = 0;
                for( size_t k = 0; k < m.size(); ++k )
                {
                    if( m[k][i].value == 1 )
                        ++onesOnFst;

                    if( m[k][j].value == 1 )
                        ++onesInSnd;
                }

                if( onesOnFst >= onesInSnd )
                {
                    answ.push_back( i );
                    break;
                }
            }
        }
    
    return answ;
}
