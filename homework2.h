#pragma once
#include <iostream>
#include <fstream>
#include "Vector.h"
#include "Matrix.h"
#include "QR_decomposition.h"
#include "Gauss.h"

MatrixD matrix_gen();
int homework2();
using namespace std;

MatrixD matrix_gen() {
	// 生成矩阵
	MatrixD m(10, 10);
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			if (i == j) {
				m[i][j] = 1.52 * cos(double(i) + 1. + 1.2 * (double(j) + 1.));
			}
			else {
				m[i][j] = sin(0.5 * i + 0.2 * j + 0.7);
			}
		}
	}
	return m;
}

int homework2() {
	fstream outfile;
	outfile.open("data2.txt", ios::out | ios::trunc);
	outfile.setf(ios::scientific, ios::floatfield);
	outfile << setprecision(12);
	MatrixD m = matrix_gen();
	QRDecomposition qr_dcmp(m);
	MatrixD qsA = qr_dcmp.getA();
	qr_dcmp.double_step_decomposition();
	MatrixD A = qr_dcmp.getA();
	VectorD er = qr_dcmp.er;
	VectorD ei = qr_dcmp.ei;

	outfile << "拟上三角化后矩阵为: \n";
	file_printmat(outfile, qsA);
	outfile << "双步位移QR分解后最终矩阵为: \n";
	file_printmat(outfile, A);
	outfile << "特征值为:\n";
	file_printvec(outfile, er);
	file_printvec(outfile, ei);
	outfile << "特征值对应的特征向量为:\n";

	int n = A.ncol();
	for (int i = 0; i < 10; i++) {
		// 如果为实根则输出特征向量
		if (ei[i] == 0) {
			double real_eigenvalue = er[i];
			VectorD d(n, real_eigenvalue);
			MatrixD to_d = diag(d);
			to_d = m - to_d;

			VectorD b(n, 0.);
			gauss(to_d, b);
			to_d[n - 1][n - 1] = 1;
			b[n - 1] = 1;
			VectorD x = back_substitution(to_d, b);
			outfile << er[i] << "\n";
			file_printvec(outfile, x);
		}
	}
	return 0;
}
