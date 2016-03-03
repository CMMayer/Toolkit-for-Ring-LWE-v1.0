#include "Transformations/TransformationRealPrimePower.h"
#include "Transformations/MatrixRealGS.h"
#include "Transformations/MatrixRealSampleGauss.h"
#include "Transformations/Transformation algorithms.h"

using RLWE_Toolkit::Transformations::TransformationRealPrimePower;

typedef TransformationRealPrimePower::TransformationType TransformationType;
typedef TransformationRealPrimePower::entry_type entry_type;

TransformationRealPrimePower::TransformationRealPrimePower(int m, int p, TransformationType transformation) :
	AbstractVectorTransformation((p - 1) * (m / p)) // dimension is phi(m) = (p-1) m' = (p-1) m/p
{
	int m_prime = m / p;

	switch (transformation)
	{
	case TransformationType::REAL_GS_D:
		matrix_p_ = std::make_unique<MatrixRealGS>(p - 1, m_prime, true);
		break;
	case TransformationType::REAL_GS_U:
		matrix_p_ = std::make_unique<MatrixRealGS>(p - 1, m_prime, false);
		break;
	case TransformationType::REAL_CONVERT_GAUSSIANS:
		matrix_p_ = std::make_unique<MatrixRealSampleGauss>(p - 1);
		break;
	}
}

TransformationRealPrimePower::~TransformationRealPrimePower(){}

// Apply this matrix to @param vec and store result in x.
void TransformationRealPrimePower::applyToVector(entry_type_vec& x, entry_type_vec const& vec)const
{
	// dim_ = phi(m), matrix_p_->getDim() = (p-1) = phi(p), m' = phi(m)/phi(p)
	int m_prime = dim_ / matrix_p_->getDim();
	applySingleKroneckerDecomposedMatrix<entry_type>(x, matrix_p_.get(), vec, 1, m_prime);
}

