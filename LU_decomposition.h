#pragma once
#include "commons.h"

class LUDecomposition {
private:
	int s = NULL;// ¾ØÕó´ø¿í, Ä¬ÈÏ²»Ï¡Êè
	int n = NULL; // ¾ØÕó½×Êý
	MatrixD m;
	VectorD b;
	VectorD x;

public:
	LUDecomposition();
	VectorD solve(const MatrixD&, const VectorD&);
	void dcmp(MatrixD&, int);
	VectorD back_substitution(const MatrixD&, const VectorD&, int);
	void set_bandwidth(const int);
	~LUDecomposition();
};

LUDecomposition::LUDecomposition() {
}

VectorD LUDecomposition::solve(const MatrixD& mat, const VectorD& bb) {
	m = mat;
	n = m.nrow();
	b = bb;
	x = VectorD(n, 0.);

	if (s == NULL) s = n;
	dcmp(m, s);
	x = back_substitution(m, b, s);
	return x;
}

void LUDecomposition::dcmp(MatrixD& mat, int s = -1) {
	if (s <= 0) s = mat.ncol();
	int n = mat.ncol();

	int temp1, temp2;
	for (int k = 0; k < n; k++) {
		temp1 = min(n, k + s + 1);
		temp2 = max(0, k - s);
		for (int j = k; j < temp1; j++) {
			for (int t = temp2; t < k; t++) {
				mat[k][j] -= mat[k][t] * mat[t][j];
			}
		}

		for (int i = k + 1; i < temp1; i++) {
			for (int t = temp2; t < k; t++) {
				mat[i][k] -= mat[i][t] * mat[t][k];
			}
			mat[i][k] = mat[i][k] / mat[k][k];
		}
	}
}

VectorD LUDecomposition::back_substitution(const MatrixD& mat, const VectorD& b, int s = -1) {
	if (s <= 0) s = mat.ncol();
	int n = mat.ncol();
	int temp = 0;
	double temp_sum = 0;
	VectorD result = b;

	// Ly = b
	for (int i = 1; i < n; i++) {
		temp = max(0, i - s);
		temp_sum = 0;
		for (int t = temp; t < i; t++) {
			temp_sum -= mat[i][t] * result[t];
		}
		result[i] = result[i] + temp_sum;
	}

	// Ux = y
	result[n - 1] = result[n - 1] / mat[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--) {
		temp = min(i + s + 1, n);
		temp_sum = 0;
		for (int t = i + 1; t < temp; t++) {
			temp_sum -= mat[i][t] * result[t];
		}
		result[i] = (result[i] + temp_sum) / mat[i][i];
	}
	return result;
}

void LUDecomposition::set_bandwidth(const int bw) {
	s = bw;
}

LUDecomposition::~LUDecomposition() {
}


