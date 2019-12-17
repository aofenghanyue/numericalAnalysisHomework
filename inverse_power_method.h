#ifndef INVERSE_POWER_METHOD_H
#define INVERSE_POWER_METHOD_H

#include<iostream>
#include"common_tools.h"
#include"LU_decomposition_old.h"
#include<vector>
#include"ctime"

template <typename T> T inverse_power_method(vector<vector<T> >& mtr, int s = -1) {
	/*���ݷ�
	*sΪ����
	*/

	T min_eigenvalue = 0;

	int row = mtr.size();
	int col = mtr[0].size();

	if (s < 0) s = row;

	// �м����u_k, y_k, beta_k
	T beta_k = 0, beta_k_1 = 0;
	vector<T > u_k(col, 1.0);
	vector<T > y_k(col, 1.0);

	// ���������������
	double err_max = 1e-7;
	double epsilon = 1;
	T u_k_norm_inf = 0;
	int norm_index = 0;
	int sgn_hr = 1;

	//clock_t begin = clock();
	// LU�ֽ�
	LU_decomposition(mtr, s);
	/*clock_t end = clock();
	cout << "LU�ֽ���ʱ" << double(end - begin) << "ms" << endl;*/

	while (epsilon > err_max) {
		// ��u_k�����
		u_k_norm_inf = norm_inf(u_k, norm_index);
		sgn_hr = sgn(u_k[norm_index]);
		// ������һ������ֵ�����ڼ������
		beta_k_1 = beta_k;
		// ��һ��
		y_k = u_k / u_k_norm_inf;

		// ���A*uk = yk
		u_k = back_substitution(mtr, y_k);//�����������������������滻��ɿ�
		//int temp = 0;
		//double temp_sum = 0;

		//for (int i = 1; i < row; i++) {
		//	temp = max(0, i - s);
		//	temp_sum = 0;
		//	for (int t = temp; t < i; t++) {
		//		temp_sum -= mtr[i][t] * y_k[t];
		//	}
		//	y_k[i] = y_k[i] + temp_sum;
		//}

		//// Ux = y
		//y_k[row - 1] = y_k[row - 1] / mtr[row - 1][col - 1];
		//for (int i = row - 2; i >= 0; i--) {
		//	temp = min(i + s + 1, col);
		//	temp_sum = 0;
		//	for (int t = i + 1; t < temp; t++) {
		//		temp_sum -= mtr[i][t] * y_k[t];
		//	}
		//	y_k[i] = (y_k[i] + temp_sum) / mtr[i][i];
		//}

		//for (int i = 0; i < row; i++) {
		//	u_k[i] = y_k[i];
		//}
		// �滻��ֹ����

		beta_k = sgn_hr * u_k[norm_index];

		epsilon = fabs(beta_k - beta_k_1) / fabs(beta_k);
	}

	min_eigenvalue = 1 / beta_k;

	return min_eigenvalue;

}

#endif // !INVERSE_POWER_METHOD_H
