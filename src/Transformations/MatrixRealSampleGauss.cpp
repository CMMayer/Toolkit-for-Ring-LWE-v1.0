#include "Transformations/MatrixRealSampleGauss.h"
#include "Math util.h"
#include <algorithm>
#include <numeric>

#define REAL_TYPE_ZERO real_type(0)
#define REAL_TYPE_ONE real_type(1)

MatrixRealSampleGauss::MatrixRealSampleGauss(int dim)
	:
	AbstractVectorTransformation(dim)
{
	// matrix represented as vector of size dim^2
	matrix_ = std::make_unique<entry_type_matrix>(dim*dim);
	int p = dim + 1;
	real_type sqrt_2 = sqrt(2.0);
	complex_type complexRootOfUnity = computeRootOfUnity(p), power1 = complex_type(REAL_TYPE_ONE), power2;
	// compute C*B' according to the explanations from [May15]
	for (int i = 0; i < dim; i++)
	{
		power2 = power1;
		for (int j = 0; j < dim/2; j++)
		{
			(*matrix_)[i*dim + j] = sqrt_2 * power2.real();
			(*matrix_)[i*dim + dim - (j+1)] = sqrt_2 * power2.imag();
			power2 = power2 * power1;
		}
		power1 = power1 * complexRootOfUnity;
	}
}

MatrixRealSampleGauss::MatrixRealSampleGauss(MatrixRealSampleGauss const& mr)
	:
	AbstractVectorTransformation(mr.dim_),
	matrix_(std::make_unique<entry_type_matrix>(*mr.matrix_))
{}

MatrixRealSampleGauss::~MatrixRealSampleGauss(){}

// Apply this matrix to @param vec and store result in x.
void MatrixRealSampleGauss::applyToVector(entry_type_vec& x, entry_type_vec const& vec) const
{
	entry_type_vec temp = vec;
	// Usual matrix-vector multiplication via inner products.
	for (int i = 0; i < dim_; i++)
	{
		entry_type_matrix::iterator it = matrix_->begin() + (i*dim_);
		x[i] = std::inner_product(it, it + dim_, temp.begin(), 0);
	}
}

