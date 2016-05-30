#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <set>
#include <cmath>

#include "Matrix.h"

//#define DEBUG_MODE

Ellement_t::Ellement_t(size_t row, size_t col, unsigned char val)
{
	this->row = row;
	this->col = col;
	this->value = val;
}

///Read matrix from .csv file with delims like "," and "\n"
Matrix::Matrix(std::string filePath)
{
	std::ifstream ifS(filePath);
	if (!ifS.is_open())
		throw std::runtime_error(filePath + " is not valid file path");

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

			Ellement_t newEll(i, j, tempU);

			row.push_back(newEll);
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
        //throw std::logic_error("Empty matrix ");
        return;

    ListOfIndexes_t kernel = getKernelRows(mPrepared);
    // запоминаем ядерные строки, добавляя их в конец уже имеющегося покрытия
    for( const auto & r : kernel )
        rowsInCovering.push_back( mPrepared[r][0].row );

#ifdef DEBUG_MODE
	std::cerr << "\nBefore deleteRowsAndCoveredColomns( rowsInCovering );\n";
	std::cerr << "rowsInCovering: " << std::endl;

	for (auto &col : rowsInCovering)
		std::cerr << (col + 1) << '\t';
	std::cerr << std::endl;
	printMatrix(mPrepared, std::cerr);
#endif //DEBUG_MODE

    deleteRowsAndCoveredColumns(kernel, mPrepared);

#ifdef DEBUG_MODE
	std::cerr << "\nAfter deleteRowsAndCoveredColomns(rowsInCovering);\n";
	printMatrix(mPrepared, std::cerr);
#endif //DEBUG_MODE
}

void Matrix::deleteRowsAndCoveredColumns(
	Matrix::ListOfIndexes_t& r,
	Matrix_t &m
)/// удаляет строки из матрицы, строки, индексы которых указаны в списке индексов р,
 /// а также столбцы, которые этими строками покрывались
{
	size_t prepSize = m.size();
	std::sort(r.begin(), r.end());

	Matrix::ListOfIndexes_t deleteRows = r;
	Matrix::ListOfIndexes_t deleteColumns;

	for (const auto &row : r)// проходимся по индексам строк
	{
		for (size_t j = 0; j < m[row].size(); ++j)//идем по строчке
			if (m[row][j].value == 1)
				deleteColumns.push_back(j);// запоминаем столбцы которые хотим удалить
	}

	for (int i = deleteRows.size() - 1; i >= 0; --i)
		std::swap(m[deleteRows[i]], m[prepSize - i - 1]);// удаляемые строчки в конец 

	m.erase(m.begin() + (prepSize - deleteRows.size()),
		m.end()); //удаляем строчик ,которые мы переместили вниз

	std::set<ListOfIndexes_t::value_type> uniqueColomns(deleteColumns.begin(), deleteColumns.end());
	deleteColumns = ListOfIndexes_t(uniqueColomns.begin(),
		uniqueColomns.end());//Отсортировали и удалили одинаковые индексы.

#ifdef DEBUG_MODE
    if( !mPrepared.empty() )
    {
        std::cerr << "deleteColumns: " << std::endl;

        for (auto &col : deleteColumns)
            std::cerr << (mPrepared[0][col].col + 1) << '\t';
    }
#endif

	this->deleteColumns(deleteColumns, mPrepared);
}

void Matrix::printMatrix(const Matrix_t & m, std::ostream & oStream)
{
	if (m.empty() || m.at(0).empty())
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

void Matrix::deleteColumns(Matrix::ListOfIndexes_t &c, Matrix_t &m)
{
	size_t i = 0;
	std::sort(c.begin(), c.end());

	for (const auto &col : c)// идем по каждому из индесов удаляемых столбцов
	{
		/// удаляем в каждой строчке заданный столбец.
		/// Но при удалении они будут сдвигаться, поэтому приходится делать -i,
		/// где i считает количество удалённых столбцов.
		for (auto &row : m)
			row.erase(row.begin() + (col - i));
		++i;
	}
}

Matrix::ListOfIndexes_t Matrix::getKernelRows(const Matrix_t& m)
{
	Matrix::ListOfIndexes_t answ;
	for (size_t j = 0; j < m.at(0).size(); ++j)// идем по столбцу
	{
		size_t onesCount = 0;// счетчик 1 в столбце 
		size_t curPreparedRow = -1;
        for (size_t i = 0; i < m.size(); ++i)// по самому столбцу
		{
			if (m[i][j].value == 1)
			{
				++onesCount;
				curPreparedRow = i;
			}
		}

		if (onesCount == 1)
			answ.push_back(curPreparedRow);
	}

    auto ker = std::set<ListOfIndexes_t::value_type>( answ.begin(), answ.end() );
    answ = ListOfIndexes_t( ker.begin(), ker.end() );
	return answ;// запоминаем индекс ядерной строки
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
	if (mPrepared.empty())
		return;

	Matrix::ListOfIndexes_t reducingColomns = getReducingColumns(mPrepared);

#ifdef DEBUG_MODE
	std::cerr << "\n Reducing colomns:\n";
	for (const auto &reduce : reducingColomns)
		std::cerr << (mPrepared[0][reduce].col + 1) << '\t';
	std::cerr << std::endl;
#endif //DEBUG_MODE

	if (reducingColomns.empty())
		return;

	this->deleteColumns(reducingColomns, mPrepared);

#ifdef DEBUG_MODE
	std::cerr << "\nAfter deleting reducing colomns:\n";
	printMatrix(mPrepared, std::cerr);
#endif //DEBUG_MODE

}

Matrix::ListOfIndexes_t Matrix::getReducingColumns(const Matrix_t &m)
{
	Matrix::ListOfIndexes_t answ;

	for (size_t i = 0; i < m.at(0).size(); ++i)// сравнивам итый и джитый столбцец, и не равно джи
		for (size_t j = 0; j < m.at(0).size(); ++j)
		{
			if (i == j)
				continue;

			bool covering = true;
			for (size_t k = 0; k < m.size(); ++k)
				if (m[k][j].value == 1 && m[k][i].value != 1)
					//проверяем,чтобы не было случая когда в итом столбце 0, а в джитом 1
				{
					covering = false;
					break;
				}

			if (covering)
			{
				size_t onesOnFst = 0, onesInSnd = 0;// считаем единички в каждом столбце
				for (size_t k = 0; k < m.size(); ++k)
				{
					if (m[k][i].value == 1)
						++onesOnFst;

					if (m[k][j].value == 1)
						++onesInSnd;
				}

				if (onesOnFst >= onesInSnd)
				{
					answ.push_back(i);
					break;
				}
			}
		}

	return answ;
}


Matrix::Matrix_t Matrix::reduceAll()
{
	if (!mPrepared.empty() && mPrepared.at(0).empty())
		mPrepared.clear();

    size_t curSize = mPrepared.size() + (mPrepared.empty() ? 0 : mPrepared[0].size());
	//Максимальный размер матрицы в данный момент
    size_t lastSize = curSize + 1;

    prepare();

    while( !mPrepared.empty() )
	{
        //prepare();
		reduceAsColumns();
		reduceAsRows();
        prepare();

        if (!mPrepared.empty() && mPrepared.at(0).empty())
            mPrepared.clear();

        lastSize = curSize;
        curSize = mPrepared.size() + (mPrepared.empty() ? 0 : mPrepared[0].size());

        if( lastSize == curSize && curSize > 0 )
        {
            deleteRowWithMinimalOnesCount( mPrepared );
            prepare();
        }
	}

	return mPrepared;
}

void Matrix::deleteRows(ListOfIndexes_t &rows, Matrix_t &m)
{
	auto unique = std::set<ListOfIndexes_t::value_type>(rows.begin(), rows.end());
	rows = ListOfIndexes_t(unique.begin(), unique.end());
	for (size_t i = 0; i < rows.size(); ++i )
		m.erase(m.begin() + (rows[i] - i));
	
#ifdef DEBUG_MODE
	std::cerr << "\nAfter deleting rows:\n";
	this->printMatrix(mPrepared, std::cerr);
#endif //DEBUG_MODE
}

void Matrix::reduceAsRows()
{
	if (mPrepared.empty())
		return;

    /*Matrix::ListOfIndexes_t ker = getKernelRows(mPrepared);
	this->rowsInCovering.insert(rowsInCovering.end(), ker.begin(), ker.end());
    this->deleteRowsAndCoveredColumns(ker, mPrepared);*/

	Matrix::ListOfIndexes_t reducingRows = getReducingRows(mPrepared);

#ifdef DEBUG_MODE
	std::cerr << "\n Reducing rows:\n";
	for (const auto &reduce : reducingRows)
		std::cerr << (mPrepared[reduce][0].row + 1) << '\t';
	std::cerr << std::endl;
#endif //DEBUG_MODE

	if (reducingRows.empty())
		return;

	deleteRows(reducingRows, mPrepared);

#ifdef DEBUG_MODE
	std::cerr << "\nAfter deleting reducing rows:\n";
	printMatrix(mPrepared, std::cerr);
#endif //DEBUG_MODE
}

Matrix::ListOfIndexes_t Matrix::getReducingRows(const Matrix_t &m)
{
	Matrix::ListOfIndexes_t answer;
	if (m.empty() || m[0].empty())
	{
		//m.clear();
		return answer;
	}
	for (size_t i = 0; i < m.size(); ++i)
		for (size_t j = i + 1; j < m.size(); ++j)
		{
			bool equals = true;

			for( size_t k = 0; k < m.at(0).size(); ++k )
				if (m[i][k].value != m[j][k].value)
				{
					equals = false;
					break;
				}

			if (equals)
				answer.push_back(i);
		}

	for (size_t i = 0; i < m.size(); ++i)
		for (size_t j = 0; j < m.size(); ++j)
		{
			if (i == j || std::find(answer.begin(), answer.end(), i) != answer.end() )
				continue;

			bool covering = true;

			for (size_t k = 0; k < m.at(0).size(); ++k)
				if (m[i][k].value == 1 && m[j][k].value != 1)
				{
					covering = false;
					break;
				}

			if (covering)
			{
				size_t onesOnFst = 0, onesInSnd = 0;
				for (size_t k = 0; k < m.at(0).size(); ++k)
				{
					if (m[i][k].value == 1)
						++onesOnFst;

					if (m[j][k].value == 1)
						++onesInSnd;
				}

				if (onesOnFst < onesInSnd)
				{
					answer.push_back(i);
					break;
				}
			}
		}

	return answer;
}

bool Matrix::isFullCovering( const ListOfIndexes_t &rows, const Matrix_t &m )
{
    if( m.empty() )
        return true;

    for( size_t j = 0; j < m.at(0).size(); ++j )
        if( std::find( rows.begin(), rows.end(), m[j][0].col) == rows.end() )
            return false;

    return true;
}

Matrix::ListOfIndexes_t Matrix::getCurCovering()
{
    std::sort( rowsInCovering.begin(), rowsInCovering.end() );
    return rowsInCovering;
}

void Matrix::deleteRowWithMinimalOnesCount( Matrix_t &m)
{
    std::vector<size_t> onesCnts( m.size(), 0u );

    for( size_t i = 0; i < m.size(); ++i )
         onesCnts[i] = std::count_if( m[i].begin(), m[i].end(),
                                      []( const Matrix_t::value_type::value_type & el) -> bool
                                      {
                                        return el.value == 1;
                                      });
    size_t indMin = std::min_element( onesCnts.begin(), onesCnts.end() ) - onesCnts.begin();
#ifdef DEBUG_MODE
    std::cerr << "\nDelete " << (m[indMin][0].row + 1 ) << " row.\n";
#endif //DEBUG_MODE
    m.erase( m.begin() + indMin );
}
