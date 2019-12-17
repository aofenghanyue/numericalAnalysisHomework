#pragma once
#include "Gauss.h"
#include "commons.h"

class MultidimensionNewtonRoot {
private:
	int n;
public:
	MatrixD jac_mat(void *jac(VectorD &x));
};
