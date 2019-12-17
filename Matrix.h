#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdio.h>

using namespace std;

template <class T>
class Matrix {
private:
	int row;
	int col;
	T** elements;
public:
	Matrix();
	Matrix(int n, int m);
	Matrix(int n, int m, const T& a);
	Matrix(int n, int m, const T* a);
	Matrix(const Matrix& other_matrix);
	Matrix(int n, int m, const Matrix& other_matrix);
	Matrix& operator=(const Matrix& other_matrix);
	inline T* operator[](const int i) const;
	inline int nrow() const;
	inline int ncol() const;
	inline Matrix<T> t();
	~Matrix();
};


template<class T>
Matrix<T>::Matrix() : row(0), col(0), elements(NULL) {}

template<class T>
Matrix<T>::Matrix(int n, int m) : row(n), col(m), elements(new T* [n]) {
	int matrix_size = n * m;
	elements[0] = new T[matrix_size];
	for (int i = 1; i < n; i++) {
		elements[i] = elements[i - 1] + m;
	}
}

template<class T>
Matrix<T>::Matrix(int n, int m, const T& a) : row(n), col(m), elements(n > 0 ? new T * [n] : NULL) {
	int matrix_size = n * m;
	if (elements != NULL) {
		elements[0] = new T[matrix_size];
		for (int i = 1; i < n; i++) {
			elements[i] = elements[i - 1] + m;
		}
		for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) elements[i][j] = a;
	}
}

template<class T>
Matrix<T>::Matrix(int n, int m, const T* a) : row(n), col(m), elements(new T* [n]) {
	int matrix_size = n * m;
	elements[0] = new T[matrix_size];
	for (int i = 1; i < n; i++) {
		elements[i] = elements[i - 1] + m;
	}
	for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) elements[i][j] = *a++;
}

template<class T>
Matrix<T>::Matrix(const Matrix& other_matrix) : row(other_matrix.row), col(other_matrix.col), elements(new T* [row]) {
	int matrix_size = row * col;
	elements[0] = new T[matrix_size];
	for (int i = 1; i < row; i++) {
		elements[i] = elements[i - 1] + col;
	}
	for (int i = 0; i < row; i++) for (int j = 0; j < col; j++) elements[i][j] = other_matrix[i][j];
}

template<class T>
inline Matrix<T>::Matrix(int n, int m, const Matrix& other_matrix) : row(n), col(m), elements(new T* [row]) {
	int matrix_size = row * col;
	elements[0] = new T[matrix_size];
	for (int i = 1; i < row; i++) {
		elements[i] = elements[i - 1] + col;
	}
	for (int i = 0; i < row; i++) for (int j = 0; j < col; j++) elements[i][j] = other_matrix[i][j];

}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix& other_matrix) {
	if (row != other_matrix.row || col != other_matrix.col) {
		if (elements) {
			delete[] elements[0];
			delete[] elements;
		}
		elements = new T * [other_matrix.row];
		this->row = other_matrix.row;
		this->col = other_matrix.col;
		int matrix_size = row * col;
		elements[0] = new T[matrix_size];
		for (int i = 1; i < row; i++) {
			elements[i] = elements[i - 1] + col;
		}
	}
	for (int i = 0; i < row; i++) for (int j = 0; j < col; j++) elements[i][j] = other_matrix[i][j];
	return *this;
}

template<class T>
inline T* Matrix<T>::operator[](const int i) const {
	return elements[i];
}

template<class T>
inline int Matrix<T>::nrow() const {
	return row;
}

template<class T>
inline int Matrix<T>::ncol() const {
	return col;
}

template<class T>
inline Matrix<T> Matrix<T>::t() {
	Matrix<T> A_t(col, row);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			A_t[j][i] = elements[i][j];
		}
	}
	return A_t;
}

template<class T>
Matrix<T>::~Matrix() {
	if (elements) {
		delete[] elements[0];
		delete[] elements;
	}
}

typedef Matrix<double> MatrixD;

#endif