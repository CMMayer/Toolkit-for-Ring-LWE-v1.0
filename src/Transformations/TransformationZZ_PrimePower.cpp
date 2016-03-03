#include "Transformations/TransformationZZ_PrimePower.h"
#include "Transformations/Transformation algorithms.h"
#include "Transformations/MatrixZZ.h"

using RLWE_Toolkit::Transformations::TransformationZZ_PrimePower;

typedef TransformationZZ_PrimePower::TransformationType TransformationType;
typedef TransformationZZ_PrimePower::entry_type entry_type;

TransformationZZ_PrimePower::TransformationZZ_PrimePower(int m, int p, TransformationType transformation) :
	AbstractVectorTransformation((p - 1) * (m / p)) // dimension is phi(m) = (p-1) m' = (p-1) m/p
{
	typedef MatrixZZ::MatrixType mat_type;
	int m_prime = m / p;

	switch (transformation)
	{
	case TransformationType::INT_L:
		matrix_p_ = std::make_unique<MatrixZZ>(p - 1, mat_type::MATRIX_L);
		break;
	case TransformationType::INT_L_INVERSE:
		matrix_p_ = std::make_unique<MatrixZZ>(p - 1, mat_type::MATRIX_L_INVERSE);
		break;
	case TransformationType::INT_G_DECODING:
		matrix_p_ = std::make_unique<MatrixZZ>(p - 1, mat_type::MATRIX_G_DECODING);
		break;
	case TransformationType::INT_G_INVERSE_DECODING:
		matrix_p_ = std::make_unique<MatrixZZ>(p - 1, mat_type::MATRIX_G_INVERSE_DECODING);
		break;
	case TransformationType::INT_G_POWERFUL:
		matrix_p_ = std::make_unique<MatrixZZ>(p - 1, mat_type::MATRIX_G_POWERFUL);
		break;
	case TransformationType::INT_G_INVERSE_POWERFUL:
		matrix_p_ = std::make_unique<MatrixZZ>(p - 1, mat_type::MATRIX_G_INVERSE_POWERFUL);
		break;
	}
}

TransformationZZ_PrimePower::~TransformationZZ_PrimePower(){}

// Apply this matrix to @param vec and store result in x.
void TransformationZZ_PrimePower::applyToVector(entry_type_vec& x, entry_type_vec const& vec)const
{
	// dim_ = phi(m), matrix_p_->getDim() = (p-1) = phi(p), m' = phi(m)/phi(p)
	int m = dim_ / matrix_p_->getDim();
	applySingleKroneckerDecomposedMatrix<entry_type>(x, matrix_p_.get(), vec, 1, m);
}

