#include "Transformations/TransformationCompDFT.h"
#include "Transformations/Transformation algorithms.h"
#include "Math util.h"

using RLWE_Toolkit::Transformations::TransformationCompDFT;
using RLWE_Toolkit::Transformations::MatrixCompFFT;

TransformationCompDFT::TransformationCompDFT(int m, int p, bool reversed, std::shared_ptr<MatrixCompFFT const> DFT_p) :
	AbstractVectorTransformation(m),
	reversed_(reversed),
	DFT_p_(DFT_p)
{
	int m_prime = m / p;
	if (m_prime > 1){ // Transformation can be further decomposed
		// recursive call of constructor
		DFT_m_prime_ = std::make_unique<TransformationCompDFT>(m_prime, p, reversed, DFT_p);

		twiddleMatrix_ = std::make_unique<entry_type_vec>(m);

		// m-th rot of unity
		entry_type tempRootOfUnity = RLWE_Toolkit::Math_util::computeRootOfUnity(m);

		if (reversed_){ //compute adjoint twiddle matrix
			for (int i = 0; i < p; i++)
			{
				for (int j = 0; j < m_prime; j++)
				{
					(*twiddleMatrix_)[i*m_prime + j] = conj(pow(tempRootOfUnity, i*j));
				}
			}
		}
		else{ //compute twiddle matrix
			for (int i = 0; i < p; i++)
			{
				for (int j = 0; j < m_prime; j++)
				{
					(*twiddleMatrix_)[i*m_prime + j] = pow(tempRootOfUnity, i*j);
				}
			}
		}
	}
}

TransformationCompDFT::~TransformationCompDFT(){}

// Apply this matrix to @param vec and store result in x.
void TransformationCompDFT::applyToVector(entry_type_vec& x, entry_type_vec const& vec) const
{
	int p = DFT_p_->getDim();
	int m_prime = dim_ / p;
	if (m_prime <= 1){ // prime case (m = p): transformation is not decomposed. twiddleMatrix_ and DFT_p_ are null pointers.
		DFT_p_->applyToVector(x, vec);
	}
	else{ // apply matrices according to the sparse decomposition
		if (reversed_){
			applySingleKroneckerDecomposedMatrix<entry_type>(x, DFT_m_prime_.get(), vec, p, 1);
			std::transform(twiddleMatrix_->cbegin(), twiddleMatrix_->cend(), x.cbegin(), x.begin(), std::multiplies<entry_type>());
			applySingleKroneckerDecomposedMatrix<entry_type>(x, DFT_p_.get(), x, 1, m_prime);
		}
		else{
			applySingleKroneckerDecomposedMatrix<entry_type>(x, DFT_p_.get(), vec, 1, m_prime);
			std::transform(twiddleMatrix_->cbegin(), twiddleMatrix_->cend(), x.cbegin(), x.begin(), std::multiplies<entry_type>());
			applySingleKroneckerDecomposedMatrix<entry_type>(x, DFT_m_prime_.get(), x, p, 1);
		}
	}
}

