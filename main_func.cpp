#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "common_tools.h"
#include "power_method.h"
#include "inverse_power_method.h"
#include "LU_decomposition.h"
#include <fstream>
#include <ctime>

using namespace std;
typedef vector<vector<double> > matrix;
vector<vector<double> > creat_matrix(int, int);

int main() {
	// 打开文件
	ofstream out_file;
	out_file.setf(ios::scientific, ios::floatfield);
	out_file << setprecision(6);
	out_file.open("data.txt", ios::out | ios::trunc);

	// 定义生成矩阵的形状
	int row = 1000, col = 1000;
	int s = 3;
	matrix mat = creat_matrix(row, col);

	// 生成矩阵个数
	int matrix_num = 81;

	// 生成列表存储条件数及最大特征值
	vector<double> cond_list(matrix_num, 0);
	vector<double> real_max_eigenvalue_list(matrix_num, 0);

	for (int k = 0; k < matrix_num; k++) {

		double t_k = k * 0.05;
		out_file << "tk = " << t_k << "时：" << endl;

		clock_t begin = clock();
		//生成矩阵
		for (int i = 0; i < row; i++) {
			for (int j = max(i - s, 0); j < min(col, i + s + 1); j++) {
				if (i == j) {
					mat[i][j] = 100 * (double(i) + 1.1) / (double(i) + 1.0 + t_k);
				}
				else {
					mat[i][j] = -(double(i) + double(j) + 2.0) / (1.0 + t_k);
				}
			}
		}

		clock_t end = clock();
		out_file <<  "生成矩阵" << "用时" << double(end - begin) << "ms" << endl;

		begin = clock();
		// 计算按模最大特征值
		double max_eigenvalue = power_method(mat, s);

		end = clock();
		out_file <<  "按模最大特征值" << "用时" << double(end - begin) << "ms" << endl;


		begin = clock();
		// 计算最大特征值
		double real_max_eigenvalue = 0;
		if (max_eigenvalue >= 0) {
			real_max_eigenvalue = max_eigenvalue;
		}
		else {
			matrix_translation(mat, -max_eigenvalue);
			real_max_eigenvalue = power_method(mat, s) + max_eigenvalue;
			matrix_translation(mat, max_eigenvalue);
		}

		end = clock();
		out_file <<  "最大特征值" << "用时" << double(end - begin) << "ms" << endl;

		begin = clock();
		// 计算按模最小特征值
		double min_eigenvalue = inverse_power_method(mat, s);
		end = clock();
		out_file <<  "最小特征值" << "用时" << double(end - begin) << "ms" << endl;
		// 条件数
		double cond = fabs(max_eigenvalue / min_eigenvalue);

		// 存储结果
		out_file << "cond = " << cond << endl;
		out_file << "最大特征值为：" << real_max_eigenvalue << endl;
		cond_list[k] = cond;
		real_max_eigenvalue_list[k] = real_max_eigenvalue;
		out_file << endl;
	}

	//使条件数最大或最小的t_k
	int index = 0;
	double cond_max = get_max_from_vector(cond_list, index);
	double cond_max_t_k = index * 0.05;
	double cond_min = get_min_from_vector(cond_list, index);
	double cond_min_t_k = index * 0.05;

	out_file << "当tk = " << cond_min_t_k << "时，条件数最小，为 " << cond_min << endl;
	out_file << "当tk = " << cond_max_t_k << "时，条件数最大，为 " << cond_max << endl;
	out_file << "当t = 0~80 时，最大特征值分别为: "<< endl;
	for (int i = 0; i < matrix_num; i++) {
		out_file << real_max_eigenvalue_list[i] << '\t';
		if (i % 3 == 2) {
			out_file << endl;
		}
	}

	out_file.close();


	return 0;
}

matrix creat_matrix(int row, int col) {
	/*初始化生成row*col的0矩阵*/
	matrix p(row, vector<double>(col, 0));
	return p;
}