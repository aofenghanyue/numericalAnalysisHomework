#pragma once
#ifndef COMMON_H
#define COMMON_H

#include <tuple>
#include <cmath>
#include "Vector.h"
#include "Matrix.h"
#include <fstream>
#include <functional>
#include <algorithm>

const double eps = 1e-12;

//template <typename T> inline int sgn(T val) {
//	/*符号函数，正值传出1，负值传出-1.*/
//	return (T(0) < val) - (val < T(0));
//}

double combination(int all, int m) {
	// 求组合数C_{all}^m
	if (m > (all / 2.)) m = all - m;
	double result = 1.;
	for (int i = 0; i < m; i++) {
		result *= (all - i) / (i + 1);
	}
	return result;
}

template <typename T> T norm_inf(Vector<T>& vec, int& index) {
	/*求向量无穷范数
	*输入向量，索引指针
	*输出无穷范数，其索引赋值给指针
	*/
	int len = vec.size();
	index = 0;

	if (len == 1) {
		return vec[0];
	}

	T max = vec[0];

	for (int i = 1; i < len; i++) {
		if (fabs(vec[i]) > max) {
			max = fabs(vec[i]);
			index = i;
		}
	}
	return max;
}

double mean(const VectorD& v) {
	int n = v.size();
	double mean = 0.;
	for (int i = 0; i < n; i++) {
		mean += v[i];
	}
	mean /= n;
	return mean;
}

double variance(const VectorD& v) {
	int n = v.size();
	double mu = mean(v);
	double var = 0.;
	for (int i = 0; i < n; i++) {
		var += pow(v[i] - mu, 2);
	}
	var /= n;
	return var;
}

VectorD solve_matrix_root_2d(double a11, double a12, double a21, double a22) {
	//传入二维矩阵的四个数
	// 创建数组存储根(a1,b1i,a2,b2i)
	VectorD root(4, 0.);
	double a = 1, b, c, delta;
	b = -a11 - a22;
	c = a11 * a22 - a12 * a21;
	delta = b * b - 4 * a * c;
	if (delta >= 0) {
		root[0] = (-b - sqrt(delta)) / (2 * a);
		root[2] = (-b + sqrt(delta)) / (2 * a);
	}
	else {
		root[0] = -b / (2 * a);
		root[1] = sqrt(-delta) / (2 * a);
		root[2] = -b / (2 * a);
		root[3] = -sqrt(-delta) / (2 * a);
	}
	return root;
}


template <typename T>
inline Vector<T> operator*(Matrix<T>& m, Vector<T>& v) {
	int len = m.nrow();
	int d = m.ncol();
	Vector <T> result(len);
	for (int i = 0; i < len; i++) {
		T sum = 0;
		for (int k = 0; k < d; k++) {
			sum += v[k] * m[i][k];
		}
		result[i] = sum;
	}
	return result;
}

template <typename T>
inline Vector<T> operator*(T b, Vector<T>& v) {
	int len = v.size();
	Vector <T> result(len);
	for (int i = 0; i < len; i++) {
		result[i] = v[i] * b;
	}
	return result;
}

template <typename T>
inline Vector<T> operator-(T b, Vector<T>& v) {
	int len = v.size();
	Vector <T> result(len);
	for (int i = 0; i < len; i++) {
		result[i] = b - v[i];
	}
	return result;
}

template <typename T>
inline Vector<T> operator-(Vector<T>& v, T b) {
	int len = v.size();
	Vector <T> result(len);
	for (int i = 0; i < len; i++) {
		result[i] = v[i] - b;
	}
	return result;
}

template <typename T>
inline Vector<T> operator*(Vector<T>& v, Matrix<T>& m) {
	int length = m.ncol();
	int d = m.nrow();
	Vector <T> result(length);
	for (int j = 0; j < length; j++) {
		T sum = 0;
		for (int k = 0; k < d; k++) {
			sum += v[k] * m[k][j];
		}
		result[j] = sum;
	}
	return result;
}

template<typename T>
inline Matrix<T> operator*(Vector<T>& v, Vector<T>& v2) {
	int length1 = v.size();
	int length2 = v2.size();
	Matrix <T> result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result[i][j] = v[i] * v2[j];
		}
	}
	return result;
}

template <typename T>
inline T dot(Vector<T>& v1, Vector<T>& v2) {
	int len = v1.size();
	T result = 0;
	for (int i = 0; i < len; i++) {
		result += v1[i] * v2[i];
	}
	return result;
}


template <typename T>
inline Vector<T> operator-(Vector<T>& v1, Vector<T>& v2) {
	int len = v1.size();
	Vector <T> result(len);
	for (int i = 0; i < len; i++) {
		result[i] = v1[i] - v2[i];
	}
	return result;
}

template <typename T>
inline Vector<T> operator+(Vector<T>& v1, Vector<T>& v2) {
	int len = v1.size();
	Vector <T> result(len);
	for (int i = 0; i < len; i++) {
		result[i] = v1[i] + v2[i];
	}
	return result;
}

template<typename T>
inline Matrix<T> operator-(Matrix<T>& m1, Matrix<T>& m) {
	int row = m1.nrow();
	int col = m1.ncol();
	Matrix <T> result(row, col);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			result[i][j] = m1[i][j] - m[i][j];
			if (fabs(result[i][j]) < eps) result[i][j] = 0;
		}
	}
	return result;
}

template<typename T>
inline Matrix<T> operator-(Matrix<T>& m1, const T& m) {
	int row = m1.nrow();
	int col = m1.ncol();
	Matrix <T> result(row, col);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			result[i][j] = m1[i][j] - m;
			if (fabs(result[i][j]) < eps) result[i][j] = 0;
		}
	}
	return result;
}

template<typename T>
inline Matrix<T> operator*(Matrix<T>& m1, Matrix<T>& m2) {
	int row = m1.nrow();
	int col = m2.ncol();
	int k_max = m1.ncol();
	Matrix <T> result(row, col);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			T sum = 0;
			for (int k = 0; k < k_max; k++) {
				sum += m1[i][k] * m2[k][j];
			}
			result[i][j] = sum;
		}
	}
	return result;
}

MatrixD diag(VectorD& v) {
	int n = v.size();
	MatrixD result(n, n, 0.);
	for (int i = 0; i < n; i++) {
		result[i][i] = v[i];
	}
	return result;
}

MatrixD eye(int n) {
	MatrixD result(n, n, 0.);
	for (int i = 0; i < n; i++) {
		result[i][i] = 1;
	}
	return result;
}

template<typename T>
void printvec(Vector<T>& v, string s = "") {
	if (s != "") cout << s << ":\n";
	for (int i = 0; i < v.size(); i++) {
		cout << v[i] << "\t";
	}
	cout << "\n";
}

template<typename T>
void printmat(Matrix<T>& m, string s = "") {
	if (s != "") cout << s << ":\n";
	for (int i = 0; i < m.nrow(); i++) {
		for (int j = 0; j < m.ncol(); j++) {
			cout << m[i][j] << "\t";
		}
		cout << "\n";
	}
}

template<typename T>
void print2matlab(Matrix<T>& m) {
	cout << "[";
	for (int i = 0; i < m.nrow(); i++) {
		for (int j = 0; j < m.ncol(); j++) {
			cout << m[i][j] << ",";
		}
		cout << "\b;";
		cout << "\n";
	}
	cout << "\b]\n";
}

template<typename T>
void file_printvec(fstream& outfile, Vector<T>& v) {
	int n = v.size();
	for (int i = 0; i < n; i++) {
		outfile << v[i] << "\t";
	}
	outfile << "\n";
}

template<typename T>
void file_printmat(fstream& outfile, Matrix<T>& m) {
	for (int i = 0; i < m.nrow(); i++) {
		for (int j = 0; j < m.ncol(); j++) {
			outfile << m[i][j] << "\t";
		}
		outfile << "\n";
	}
}

#endif
