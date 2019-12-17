#pragma once
#include "commons.h"
#include "LU_decomposition.h"
#include "Gauss.h"

class BivariatePolynomialFitting {
private:
	int N = 1; //x多项式最高阶数+1;
	int n = 0;
	int M = 1; //y多项式最高阶数+1;
	int m = 0;
	string method = "lu";
	double total_error = 1; //总误差
	MatrixD C; //系数矩阵c=(B_T*B)^(-1) (B_T*U <-A|||D_T->* G) (G_T*G)^(-1)
	VectorD x_nodes, y_nodes; //节点 n x 1, m x 1
	MatrixD Z; //节点值z_ij = z(x_i, y_j) (n, m)
	MatrixD B; //x的基函数在各节点的值[phi0, ...] (n, N)
	MatrixD G; //y的基函数在各节点的值[psi0, ...] (m, M)
public:
	BivariatePolynomialFitting(int, int);
	void fit(const VectorD&, const VectorD&, const MatrixD&);
	void fit_coef_lu();
	void fit_coef_gauss();
	double predict(double, double);
	double fit_error();
	double get_error();
	MatrixD get_coef();
	void set_method(string);
	MatrixD tran_coef(const MatrixD&, double mu_x, double sigma_x, double mu_y, double sigma_y);
	~BivariatePolynomialFitting();
};

BivariatePolynomialFitting::BivariatePolynomialFitting(int n, int m) : N(n), M(m) {
}

void BivariatePolynomialFitting::fit(const VectorD& x_data, const VectorD& y_data, const MatrixD& z_data) {
	/*调用此函数进行拟合
	* 输入：
	* x_data, y_data: x, y 节点
	* z_data: z_ij = z(x_i, y_j)
	*/
	// 节点个数
	n = x_data.size();
	m = y_data.size();

	// 初始化各符号
	x_nodes = x_data;
	y_nodes = y_data;
	Z = z_data;

	// 取phi_i(x)=x^i，计算B矩阵
	B = MatrixD(n, N);
	for (int ii = 0; ii < n; ii++) {
		for (int jj = 0; jj < N; jj++) {
			B[ii][jj] = pow(x_nodes[ii], jj);
		}
	}

	// 取psi_j(y)=y^j，计算G矩阵
	G = MatrixD(m, M);
	for (int ii = 0; ii < m; ii++) {
		for (int jj = 0; jj < M; jj++) {
			G[ii][jj] = pow(y_nodes[ii], jj);
		}
	}

	// 根据所选模式进行拟合
	method == "lu" ? fit_coef_lu() : fit_coef_gauss();

	// 记录误差
	total_error = fit_error();
}

void BivariatePolynomialFitting::fit_coef_lu() {
	MatrixD B_T = B.t();
	MatrixD G_T = G.t();
	MatrixD Z_T = Z.t();

	LUDecomposition lu;
	// 利用LU分解求A=(B_T*B)^(-1) B_T*U, 即(B_T*B)*A = B_T*U
	MatrixD A(N, m);
	MatrixD coef = B_T * B;
	lu.dcmp(coef);

	for (int ii = 0; ii < m; ii++) {
		VectorD b1(n, Z_T[ii]);
		b1 = B_T * b1;
		VectorD x = lu.back_substitution(coef, b1);
		for (int jj = 0; jj < N; jj++) {
			A[jj][ii] = x[jj];
		}
	}

	// 求D_T=G*(G_T*G)^(-1), 即G_T*G*D = G_T
	MatrixD D(M, m);
	MatrixD D_T;
	coef = G_T * G;
	lu.dcmp(coef);
	for (int ii = 0; ii < m; ii++) {
		VectorD b2(M, G[ii]);
		VectorD x = lu.back_substitution(coef, b2);
		for (int jj = 0; jj < M; jj++) {
			D[jj][ii] = x[jj];
		}
	}
	D_T = D.t();

	// 求C
	C = A * D_T;
}

void BivariatePolynomialFitting::fit_coef_gauss() {
	MatrixD B_T = B.t();
	MatrixD G_T = G.t();
	MatrixD Z_T = Z.t();

	// 利用高斯消去法求A,(B_T*B)*A = B_T*U
	MatrixD A(N, m);
	MatrixD coef = B_T * B;

	for (int ii = 0; ii < m; ii++) {
		MatrixD tempA = coef;
		VectorD b1(n, Z_T[ii]);
		b1 = B_T * b1;
		VectorD x = gauss_solve(tempA, b1);
		for (int jj = 0; jj < N; jj++) {
			A[jj][ii] = x[jj];
		}
	}

	// 求D_T=G*(G_T*G)^(-1), 即G_T*G*D = G_T
	MatrixD D(M, m);
	MatrixD D_T;
	coef = G_T * G;
	for (int ii = 0; ii < m; ii++) {
		MatrixD tempA = coef;
		VectorD b2(M, G[ii]);
		VectorD x = gauss_solve(tempA, b2);
		for (int jj = 0; jj < M; jj++) {
			D[jj][ii] = x[jj];
		}
	}
	D_T = D.t();

	// 求C
	C = A * D_T;
}

inline double BivariatePolynomialFitting::predict(double x, double y) {
	/*用拟合后的多项式预测在x,y点的函数值*/
	double result = 0.;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			result += C[i][j] * pow(x, i) * pow(y, j);
		}
	}
	return result;
}

double BivariatePolynomialFitting::fit_error() {
	/*求拟合误差*/
	int n = x_nodes.size();
	int m = y_nodes.size();
	long double err = 0.;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			err += powl(Z[i][j] - predict(x_nodes[i], y_nodes[j]), 2);
		}
	}
	return err;
}

inline double BivariatePolynomialFitting::get_error() {
	return total_error;
}

MatrixD BivariatePolynomialFitting::get_coef() {
	return C;
}

void BivariatePolynomialFitting::set_method(string mtd) {
	if (mtd == "gauss" || mtd == "lu") method = mtd;
}

MatrixD BivariatePolynomialFitting::tran_coef(const MatrixD& coef, double mu_x, double sigma_x, double mu_y, double sigma_y) {
	/*将z=c(x_new)^r*(y_new)^s展开，得到x^r*y^s系数*/
	MatrixD temp = coef;
	int row = coef.nrow();
	int col = coef.ncol();
	int max_order = max(row, col); // 多项式最高阶次数+1(包含0)
	MatrixD result(row, col, 0.);
	VectorD newton_coef(max_order); // 牛顿二项式展开系数
	for (int i = 0; i < max_order; i++) {
		newton_coef[i] = combination(max_order - 1, i);
	}
	// 系数除以sigma^...
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			temp[i][j] = temp[i][j] / pow(sigma_x, i) / pow(sigma_y, j);
		}
	}

	for (int r = 0; r < row; r++) {
		for (int s = 0; s < col; s++) {
			for (int i = 0; i <= r; i++) {
				for (int j = 0; j <= s; j++) {
					// Crs(x-mu_x)^r*(y-mu_y)^s 对x^i*y^j系数的影响
					result[i][j] += temp[r][s] * newton_coef[i] * pow(-mu_x, r-i) * newton_coef[j] * pow(-mu_y, s-j);
				}
			}
		}
	}
	return result;
}

BivariatePolynomialFitting::~BivariatePolynomialFitting() {
}

