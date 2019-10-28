
#ifndef TOOL_H
#define TOOL_H

#include <vector>
#include <algorithm>
using namespace std;

template <typename T> inline int sgn(T val) {
	/*符号函数，正值传出1，负值传出-1.*/
	return (T(0) < val) - (val < T(0));
}

template <typename T> T norm_inf(vector<T>& vec, int& index) {
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

template <typename T> T get_max_from_vector(vector<T>& vec, int& index) {
	/*求向量中最大值及其索引*/
	int len = vec.size();
	index = 0;

	if (len == 1) {
		return vec[0];
	}

	T max = vec[0];

	for (int i = 1; i < len; i++) {
		if (vec[i] > max) {
			max = vec[i];
			index = i;
		}
	}
	return max;
}

template <typename T> T get_min_from_vector(vector<T>& vec, int& index) {
	/*求向量中最小值及其索引*/
	int len = vec.size();
	index = 0;

	if (len == 1) {
		return vec[0];
	}

	T min = vec[0];

	for (int i = 1; i < len; i++) {
		if (vec[i] < min) {
			min = vec[i];
			index = i;
		}
	}
	return min;
}

template <typename T>
vector<T> operator/(vector<T>& v, T b) {
	int len = v.size();
	vector <T> result(len);
	for (int i = 0; i < len; i++) {
		result[i] = v[i] / b;
	}
	return result;
}

template <typename T>
inline vector<T> dot(vector<vector<T> >& mtr, vector<T>& vec, int row, int col, int s = -1) {
	/*带状矩阵乘向量*/
	if (s < 0) s = row;
	vector <T> result(row);
	for (int i = 0; i < row; i++) {
		T resulti = 0;
		for (int j = max(0, i - s); j < min(col, i + s + 1); j++) {
			resulti += mtr[i][j] * vec[j];
		}
		result[i] = resulti;
	}
	return result;
}

template <typename T, typename T1>
inline void matrix_translation(vector<vector<T> >& mtr, T1 a) {
	int shape = mtr.size();
	for (int i = 0; i < shape; i++) {
		mtr[i][i] += a;
	}
}
#endif // TOOL_H


