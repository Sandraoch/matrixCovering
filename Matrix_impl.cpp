#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <set>
#include <cmath>

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
		throw std::logic_error("Empty matrix ");

	auto kernel = getKernelRows(mPrepared);
	rowsInCovering.insert(rowsInCovering.end(),
		kernel.begin(),
		kernel.end()
	);// запоминаем ядерные строки, добавляя их в конец уже имеющегося покрытия 

#ifdef DEBUG_MODE
	std::cerr << "\nBefore deleteRowsAndCoveredColomns( rowsInCovering );\n";
	std::cerr << "rowsInCovering: " << std::endl;

	for (auto &col : rowsInCovering)
		std::cerr << (col + 1) << '\t';
	std::cerr << std::endl;
	printMatrix(mPrepared, std::cerr);
#endif //DEBUG_MODE

	deleteRowsAndCoveredColumns(rowsInCovering, mPrepared);

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
	std::cerr << "deleteColumns: " << std::endl;

	for (auto &col : deleteColumns)
		std::cerr << (col + 1) << '\t';
#endif

	this->deleteColumns(deleteColumns, mPrepared);
}

void Matrix::printMatrix(const Matrix_t & m, std::ostream & oStream)
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
		for (size_t i = 0; i < m.size(); i++)// по самому столбцу
		{
			if (m[i][j].value == 1 )
			{
				++onesCount;
				curPreparedRow = i;
			}
		}

		if (onesCount == 1)
			answ.push_back(curPreparedRow);
	}

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
	size_t curSize = mPrepared.size() + 1;//количество строк в данный момент

	this->prepare();

	while (curSize > mPrepared.size())
	{
		curSize = mPrepared.size();

		reduceAsColumns();
		reduceAsRows();
	}

	return mPrepared;
}

void Matrix::reduceAsRows()
{
	Matrix::ListOfIndexes_t reducingRows = getReducingRows(mPrepared);

#ifdef DEBUG_MODE
	std::cerr << "\n Reducing rows:\n";
	for (const auto &reduce : reducingRows)
		std::cerr << (mPrepared[reduce][0].row + 1) << '\t';
	std::cerr << std::endl;
#endif //DEBUG_MODE

	if (reducingRows.empty())
		return;

	this->rowsInCovering.insert(rowsInCovering.end(), reducingRows.begin(), reducingRows.end());
	this->deleteRowsAndCoveredColumns(reducingRows, mPrepared);

#ifdef DEBUG_MODE
	std::cerr << "\nAfter deleting reducing rows:\n";
	printMatrix(mPrepared, std::cerr);
#endif //DEBUG_MODE
}

Matrix::ListOfIndexes_t Matrix::getReducingRows(const Matrix_t &m)
{
	Matrix::ListOfIndexes_t answer;
	for (size_t i = 0; i < m.size(); ++i)
		for (size_t j = 0; j < m.size(); ++j)
		{
			if (i == j)
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

				if (onesOnFst <= onesInSnd)
				{
					answer.push_back(i);
					break;
				}
			}
		}

	return answer;
}
