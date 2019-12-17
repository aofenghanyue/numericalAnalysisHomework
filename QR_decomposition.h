#pragma once
#ifndef QR_DECOMPOSITION_H
#define QR_DECOMPOSITION_H
#include "commons.h"

using namespace std;

class QRDecomposition {
	/*QR�ֽ���*/
private:
	int n; // �������
	MatrixD Q, A; // һ��QR�ֽ��м����
	MatrixD gen_mk(MatrixD& m);
public:
	int L = 10000; // ����������
	VectorD er, ei; // eigenvalues_real, eigenvalues_imaginary, ����ֵ;
	QRDecomposition(MatrixD& a, bool initial = true);
	void quasi_triangulation();
	void double_step_decomposition();
	tuple<VectorD, double> householder_trans(int ii, int jj);
	void one_step_HAH(VectorD& ur, double& hr);
	void one_step_HA(VectorD& ur, double& hr);
	MatrixD getA();
	~QRDecomposition();
};

inline MatrixD QRDecomposition::gen_mk(MatrixD& m) {
	MatrixD msquare = m * m;
	int n = m.nrow();
	MatrixD result(n, n);
	double s = (m[n - 2][n - 2] + m[n - 1][n - 1]);
	double det = (m[n - 2][n - 2] * m[n - 1][n - 1] - m[n - 1][n - 2] * m[n - 2][n - 1]);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				result[i][j] = msquare[i][j] - s * m[i][j] + det;
			}
			else {
				result[i][j] = msquare[i][j] - s * m[i][j];
			}
		}
	}
	return result;
}


QRDecomposition::QRDecomposition(MatrixD& a, bool initial) : n(a.nrow()), A(a) {
	if (initial) quasi_triangulation();
}


void QRDecomposition::quasi_triangulation() {
	/*�������ǻ�*/
	double hr;
	VectorD ur;
	for (int r = 0; r < n - 2; r++) {
		// �жϸ����Ƿ�ȫΪ0
		for (int i = r + 2; i < n; i++) {
			if (A[i][r] != 0) {
				goto label;
			}
		}
		continue;
	label:
		// �����r����hr, ur
		tie(ur, hr) = householder_trans(r + 1, r);
		// ����Hr*A*Hr
		one_step_HAH(ur, hr);
	}
}

void QRDecomposition::double_step_decomposition() {
	er = VectorD(n);
	ei = VectorD(n);
	int m = A.nrow();
	for (int k = 0; k < L; k++) {
		if (m == 1) {
			er[0] = A[0][0];
			ei[0] = 0;
			break;
		}
		else if (m == 2) {
			VectorD root = solve_matrix_root_2d(A[0][0], A[0][1], A[1][0], A[1][1]);
			er[0] = root[0];
			er[1] = root[2];
			ei[0] = root[1];
			ei[1] = root[3];
			break;
		}
		else if (fabs(A[m - 1][m - 2]) <= eps) {
			er[m - 1] = A[m - 1][m - 1];
			ei[m - 1] = 0;
			m -= 1;
		}
		else if (fabs(A[m - 2][m - 3]) <= eps) {
			VectorD root = solve_matrix_root_2d(A[m - 2][m - 2], A[m - 2][m - 1], A[m - 1][m - 2], A[m - 1][m - 1]);
			er[m - 1] = root[0];
			er[m - 2] = root[2];
			ei[m - 1] = root[1];
			ei[m - 2] = root[3];
			m -= 2;
		}
		else {
			// ����˫��λ��QR�ֽ���м����M_k, A_k
			MatrixD A_k(m, m, A);
			MatrixD m_k = gen_mk(A_k);
			VectorD ur(m);
			double hr = 0;
			// ��m_k�����h_r
			/*cout << "\nm_k\n";
			printmat(m_k);*/
			QRDecomposition qr_mk(m_k, false);
			/*cout << "\nA_k\n";
			printmat(A_k);*/
			QRDecomposition qr_ak(A_k, false);
			for (int i = 0; i < m - 1; i++) {
				for (int kk = i + 1; i < m; i++) {
					if (A[kk][i] != 0) {
						goto label1;
					}
				}
				continue;
			label1:
				tie(ur, hr) = qr_mk.householder_trans(i, i);
				qr_ak.one_step_HAH(ur, hr);
				qr_mk.one_step_HA(ur, hr);
			}
			MatrixD A_temp = qr_ak.getA();
			for (int i = 0; i < m; i++) {
				for (int ii = 0; ii < m; ii++) {
					A[i][ii] = A_temp[i][ii];
				}
			}
		}
	}


}

tuple<VectorD, double> QRDecomposition::householder_trans(int ii, int jj) {
	/*����a��ii��(��0��ʼ),jj�п�ʼ��householder�任,������任��cr, ur, hr=1/2*||ur||_2*/
	int n = A.nrow();
	double hr, cr;
	VectorD ur(n, 0.);

	double sum = 0.;
	for (int i = ii; i < n; i++) {
		double a_ij = A[i][jj];
		sum += a_ij * a_ij;
		ur[i] = a_ij;
	}
	// sr=(0,0,...,a_rr,...,a_nr)
	double norm_sr = sqrt(sum);
	cr = -sgn(A[ii][jj]) * norm_sr;
	// ur=(0,...,a_rr-cr,...,a_nr)
	ur[ii] -= cr;
	// hr = 1/2*||ur||_2^2 = cr(cr-a_rr)
	hr = -cr * ur[ii];

	return make_tuple(ur, hr);
}

void QRDecomposition::one_step_HAH(VectorD& ur, double& hr) {
	// ����Hr*A*Hr
	MatrixD A_t = A.t();
	VectorD vr = ur / hr;
	VectorD pr = A_t * vr;
	VectorD qr = A * vr;
	double tr = dot(vr, pr);
	VectorD temp1 = tr * ur;
	VectorD wr = qr - temp1;
	for (int i = 0; i < A.nrow(); i++) {
		for (int j = 0; j < A.ncol(); j++) {
			A[i][j] -= wr[i] * ur[j] + ur[i] * pr[j];
			if (abs(A[i][j]) < eps) A[i][j] = 0;
		}
	}
}

inline void QRDecomposition::one_step_HA(VectorD& ur, double& hr) {
	// ����H*A
	MatrixD A_t = A.t();
	VectorD vr = ur / hr;
	VectorD uA = ur * A;
	MatrixD wr = vr * uA;
	for (int i = 0; i < A.nrow(); i++) {
		for (int j = 0; j < A.ncol(); j++) {
			A[i][j] -= wr[i][j];
			if (abs(A[i][j]) < eps) A[i][j] = 0;
		}
	}
}

inline MatrixD QRDecomposition::getA() {
	return A;
}

inline QRDecomposition::~QRDecomposition() {
}


#endif // !QR_DECOMPOSITION_H


