#pragma once
#include "commons.h"
#include "MultidimensionNewtonRoot.h"
#include "BivariatePolynomialFitting.h"
#include <vector>


function<MatrixD(VectorD&)> gen_jac_mat(double x, double y) {
	/*���ڴ����x,y���������Ӧ�ķ�������ſɱȾ��󣬴�������*/
	auto J = [x, y](VectorD& v) -> MatrixD {
		double m[25] = {
			-exp(-v[0]), -2 * exp(-2 * v[1]), 1, -2, x,
			-2 * exp(-2 * v[0]), -exp(-v[1]), -2, y, -1,
			x, 3, -exp(-v[2]), 0, -3,
			2, y, 1, exp(-v[3]), -4 * exp(-2 * v[4]),
			1,-2, -3 * x, -2 * exp(-2 * v[3]), -3 * exp(-v[4])
		};
		MatrixD result(5, 5, m);
		return result;
	};
	return J;
}

function<VectorD(VectorD&)> gen_fvec(double x, double y) {
	/*���ڴ����x,y���������Ӧ�ķ�������˺�������������*/
	auto f = [x, y](VectorD& v) -> VectorD {
		double vec[5] = {
			exp(-v[0]) + exp(-2 * v[1]) + v[2] - 2 * v[3] + x * v[4] - 5.3,
			exp(-2 * v[0]) + exp(-v[1]) - 2 * v[2] + y * v[3] - v[4] + 25.6,
			x * v[0] + 3 * v[1] + exp(-v[2]) - 3 * v[4] + 37.8,
			2 * v[0] + y * v[1] + v[2] - exp(-v[3]) + 2 * exp(-2 * v[4]) - 31.3,
			v[0] - 2 * v[1] - 3 * x * v[2] + exp(-2 * v[3]) + 3 * exp(-v[4]) + 42.1
		};
		VectorD result(5, vec);
		return result;
	};
	return f;
}

void homework3() {
	/*�����δ���ҵ������*/
	// �趨�����ʽ������e�ͣ�����12
	cout.setf(ios::scientific, ios::floatfield);
	cout << setprecision(12);
	// ������Ͻڵ�
	int x_num = 31;
	int y_num = 21;
	VectorD x_list(x_num, 0.);
	VectorD y_list(y_num, 0.);

	MatrixD z1(x_num, y_num, 0.);
	MatrixD z2(x_num, y_num, 0.);
	MatrixD z3(x_num, y_num, 0.);
	MatrixD z4(x_num, y_num, 0.);
	MatrixD z5(x_num, y_num, 0.);
	// ţ�ٵ�������һ��ȡ��ֵΪ����5��������Ϊ2������
	VectorD init_z(5, 2.);
	vector<MatrixD> z{ z1,z2,z3,z4,z5 };

	// ���x, y, z
	for (int i = 0; i < x_num; i++) {
		for (int j = 0; j < y_num; j++) {
			// ��x_i, y_j
			x_list[i] = 1 + 0.008 * i;
			y_list[j] = 1 + 0.008 * j;
			// ���x_i, y_j ��Ӧ���ſɱȾ�������˺�������
			function<MatrixD(VectorD&)> jac_mat = gen_jac_mat(x_list[i], y_list[j]);
			function<VectorD(VectorD&)> f_vec = gen_fvec(x_list[i], y_list[j]);
			// ţ�ٷ����z
			MultidimensionNewtonRoot newton_solve;
			newton_solve.set_precision(1e-14); // ���þ���
			newton_solve.set_jac_mat(&jac_mat); // �����ſɱȾ���
			newton_solve.set_f_vec(&f_vec); // ������������
			newton_solve.find_roots(init_z); // ��ֵȡ��һ������z����ֹ��ɢ
			VectorD z_root = newton_solve.get_roots();
			init_z = z_root;
			// �洢z
			for (int z_i = 0; z_i < 5; z_i++) {
				z[z_i][i][j] = z_root[z_i];
			}
		}
	}

	printmat(z[0], "Z1");

	// ����Ϻ�Ĳ���
	// ��׼�� x_new = (x-mu_x)/sigma_x, y_newͬ��
	double mu_x = mean(x_list);
	double mu_y = mean(y_list);
	double sig_x = sqrt(variance(x_list));
	double sig_y = sqrt(variance(y_list));
	VectorD x_list_new = (x_list - mu_x) / sig_x;
	VectorD y_list_new = (y_list - mu_y) / sig_y;
	for (int kk = 0; kk < 5; kk++) {
		cout << "z" << kk + 1 << ":\n";
		int x_poly_max = 10;
		int y_poly_max = 10;
		for (int k_t = 1; k_t < x_poly_max; k_t++) {
			for (int l_t = 1; l_t < y_poly_max; l_t++) {
				// ��С�������
				BivariatePolynomialFitting fit_model(k_t, l_t);
				fit_model.fit(x_list_new, y_list_new, z[kk]);
				double err = fit_model.get_error();

				if (err <= 1e-4) {
					// ������㾫�ȵ�ϵ������
					MatrixD coef_mat = fit_model.get_coef();
					// ���Ա�׼����x_new,y_new��ϵ������ת��Ϊ��x,y��ϵ������
					coef_mat = fit_model.tran_coef(coef_mat, mu_x, sig_x, mu_y, sig_y);
					printmat(coef_mat, "C");
					// ���k = k_t, l = l_t���㾫�ȣ���l_max��Ϊl_t
					y_poly_max = l_t;
					cout << "k = " << k_t - 1 << ", l = " << l_t - 1 << ", error=" << err << "\n";
					break;
				}

			}
		}

	}
}



