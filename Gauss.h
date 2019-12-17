#pragma once
#include "commons.h"

void gauss(MatrixD& a, VectorD& v) {
	int i, j, k, max;
	int n = a.ncol();
	double temp;
	for (k = 0; k < n - 1; k++) {
		max = k;
		for (i = k + 1; i < n; i++) {
			if (fabs(a[i][k]) > fabs(a[max][k])) {
				max = i;
			}
		}
		if (k != max) {
			for (i = k; i < n; i++) {
				temp = a[k][i];
				a[k][i] = a[max][i];
				a[max][i] = temp;
			}
			temp = v[max];
			v[max] = v[k];
			v[k] = temp;
		}

		for (i = k + 1; i < n; i++) {
			temp = a[i][k] / a[k][k];
			for (j = k; j < n; j++) {
				a[i][j] -= temp * a[k][j];
			}
			v[i] -= v[k] * temp;
		}
	}
}

VectorD back_substitution(MatrixD& a, const VectorD &b) {
	int n = a.nrow();
	VectorD x(n, 0.);
	int j, k;
	double temp;
	x[n - 1] = b[n - 1] / a[n - 1][n - 1];
	for (k = n - 2; k >= 0; k--) {
		temp = b[k];
		for (j = k + 1; j < n; j++) {
			temp -= a[k][j] * x[j];
		}
		x[k] = temp / a[k][k];
	}
	return x;
}

VectorD gauss_solve(MatrixD& A, VectorD& b) {
	gauss(A, b);
	VectorD x = back_substitution(A, b);
	return x;
}

