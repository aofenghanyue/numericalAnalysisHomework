#ifndef POWER_METHOD_H
#define POWER_METHOD_H

#include<iostream>
#include"common_tools.h"
#include<vector>

using namespace std;

template <typename T> T power_method(vector<vector<T> > &mtr, int s = -1) {
	/*幂法 
	*s为带宽
	*/
	T max_eigenvalue = 0;

	int row = mtr.size();
	int col = mtr[0].size();

	if (s < 0) s = row;

	// 中间变量u_k, y_k, beta_k
	T beta_k = 0, beta_k_1 = 0;
	vector<T > u_k(col, 1.0);
	vector<T > y_k(col, 1.0);

	// 误差，无穷范数及其索引
	double err_max = 1e-7;
	double epsilon = 1;
	T u_k_norm_inf = 0;
	int norm_index = 0;
	int sgn_hr = 1;

	T resulti = 0;
	int temp1, temp2;
	T m1, m2;

	while (epsilon > err_max) {
		// 求u_k无穷范数
		u_k_norm_inf = norm_inf(u_k, norm_index);
		sgn_hr = sgn(u_k[norm_index]);

		beta_k_1 = beta_k;

		y_k = u_k / u_k_norm_inf;

		//u_k = dot(mtr, y_k, row, col, s); //太慢，用下面的替换还是很慢
		for (int i = 0; i < row; i++) {
			resulti = 0;
			temp1 = max(0, i - s);
			temp2 = min(col, i + s + 1);
			for (int j = temp1; j < temp2; j++) {
				m1 = mtr[i][j];
				m2 = y_k[j];

				resulti += m1 * m2;
			}
			u_k[i] = resulti;
		}
		// 替换截止这里

		beta_k = sgn_hr * u_k[norm_index];

		epsilon = fabs(beta_k - beta_k_1) / fabs(beta_k);
	}

	max_eigenvalue = beta_k;

	return max_eigenvalue;

}

#endif // !POWER_METHOD_H



