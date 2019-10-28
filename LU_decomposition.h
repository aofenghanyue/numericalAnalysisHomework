#ifndef LU_DECOMPOSITION_H
#define LU_DECOMPOSITION_H

#include<iostream>
#include"common_tools.h"
#include<vector>
#include<algorithm>

template <typename T> void LU_decomposition(vector<vector<T> >& mtr, int s = -1) {
	/*LU分解*/
	int row = mtr.size();
	int col = mtr[0].size();

	if (s < 0) s = row;

	for (int k = 0; k < row; k++) {
		for (int j = k; j < min(col, k + s + 1); j++) {
			for (int t = max(0, k - s); t < k; t++) {
				mtr[k][j] -= mtr[k][t] * mtr[t][j];
			}
		}

		for (int i = k + 1; i < min(row, k + s + 1); i++) {
			for (int t = max(0, k - s); t < k; t++) {
				mtr[i][k] -= mtr[i][t] * mtr[t][k];
			}
			mtr[i][k] = mtr[i][k] / mtr[k][k];
		}
	}
}

template<typename T>
vector<T> back_substitution(const vector<vector<T> > & mtr, vector<T> &b, int s = -1) {
	/*很慢*/
	int row = mtr.size();
	int col = mtr[0].size();

	if (s < 0) s = row;
	// 回代求解
	// Ly = b
	int temp = 0;
	double temp_sum = 0;

	for (int i = 1; i < row; i++) {
		temp = max(0, i - s);
		temp_sum = 0;
		for (int t = temp; t < i; t++) {
			temp_sum -= mtr[i][t] * b[t];
		}
		b[i] = b[i] + temp_sum;
	}

	// Ux = y
	b[row - 1] = b[row - 1] / mtr[row - 1][col - 1];
	for (int i = row - 2; i >= 0; i--) {
		temp = min(i + s + 1, col);
		temp_sum = 0;
		for (int t = i + 1; t < temp; t++) {
			temp_sum -= mtr[i][t] * b[t];
		}
		b[i] = (b[i] + temp_sum) / mtr[i][i];
	}
	return b;
}

template <typename T> vector<T> LU_decomposition_solve(vector<vector<T> >& mtr, vector<T>& b, int s = -1) {
	/*LU分解求解方程*/
	int row = mtr.size();
	int col = mtr[0].size();

	if (s < 0) s = row;

	LU_decomposition(mtr, s);

	return back_substitution(mtr, b, s);
}



#endif // !LU_DECOMPOSITION_H

