#pragma once
#include "commons.h"
#include "Gauss.h"

class MultidimensionNewtonRoot {
private:
	int n;
	int iter_limit = 1000;
	double epsilon = eps;
	MatrixD jac;
	VectorD x;
	VectorD dx;
	VectorD fx;
	VectorD roots;
	function<MatrixD(VectorD&)>* jac_mat;
	function<VectorD(VectorD&)>* f_vec;
public:
	MultidimensionNewtonRoot();
	void set_f_vec(function<VectorD(VectorD&)>* f);
	void set_jac_mat(function<MatrixD(VectorD&)>* f);
	void find_roots(const VectorD& init_x);
	void set_precision(double err_max);
	void set_iter_times(int iter_time);
	VectorD get_roots();
	~MultidimensionNewtonRoot();
};

MultidimensionNewtonRoot::MultidimensionNewtonRoot() {
}

inline void MultidimensionNewtonRoot::set_f_vec(function<VectorD(VectorD&)>* f) {
	f_vec = f;
}

void MultidimensionNewtonRoot::set_jac_mat(function<MatrixD(VectorD&)>* f) {
	jac_mat = f;
}

void MultidimensionNewtonRoot::find_roots(const VectorD& init_x) {
	/*牛顿迭代法求根*/
	x = init_x;
	double err = 1.;
	while (err > epsilon) {
		// J(x)·Δx = -F(x)
		jac = (*jac_mat)(x);
		fx = (*f_vec)(x);
		fx = 0. - fx;
		dx = gauss_solve(jac, fx);
		// x_k+1 = x_k + Δx
		x = x + dx;

		int index;
		err = norm_inf(dx, index) / norm_inf(x, index);
	}
}

void MultidimensionNewtonRoot::set_precision(double err_max) {
	epsilon = err_max;
}

void MultidimensionNewtonRoot::set_iter_times(int iter_time) {
	iter_limit = iter_time;
}

VectorD MultidimensionNewtonRoot::get_roots() {
	return x;
}

MultidimensionNewtonRoot::~MultidimensionNewtonRoot() {
}
