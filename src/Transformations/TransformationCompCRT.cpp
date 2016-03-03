#include "Transformations/TransformationCompCRT.h"
#include "Transformations/TransformationCompDFT.h"
#include "Transformations/Transformation algorithms.h"
#include "Transformations/MatrixCompMult.h"

#include "Math util.h"

using RLWE_Toolkit::Transformations::TransformationCompCRT;

typedef TransformationCompCRT::TransformationType TransformationType;

TransformationCompCRT::TransformationCompCRT(int m, int p, TransformationType transformation) :
	AbstractVectorTransformation((p-1) * (m/p))
{
	typedef MatrixCompFFT::MatrixType trans_type;

	int m_prime = m / p;
	int phi_m = (p - 1) * (m / p);
	trans_type trans;

	switch (transformation){
		// switch over the indicator for the represented CRT. Construct respective prime-indexed CRT,
		// set indicator for the respective prime-indexed DFT and the reversed flag.
	case TransformationType::COMPLEX_CRT_M:
		CRT_p_ = std::make_unique<MatrixCompFFT>(trans_type::CRT_P, p);
		trans = trans_type::DFT_P;
		reversed_ = false;
		break; 

	case TransformationType::COMPLEX_CRT_M_INVERSE:
		CRT_p_ = std::make_unique<MatrixCompMult>(p, false);
		trans = trans_type::DFT_P_INV;
		reversed_ = true;
		break;

	case TransformationType::COMPLEX_CRT_M_STAR:
		CRT_p_ = std::make_unique<MatrixCompFFT>(trans_type::CRT_P_STAR, p);
		trans = trans_type::DFT_P_STAR;
		reversed_ = true;
		break;

	case TransformationType::COMPLEX_CRT_M_STAR_INVERSE:
		CRT_p_ = std::make_unique<MatrixCompMult>(p, true);
		trans = trans_type::DFT_P_STAR_INV;
		reversed_ = false;
		break;

	default:
		break;
	}
	
	if (m_prime > 1){
		// construct DFT_p
		std::shared_ptr<MatrixCompFFT> DFT_p = std::make_shared<MatrixCompFFT>(trans, p);
		// and pass it to the constructor of DFT_m'. The constructor calls itself recursively
		// until m' = p.
		DFT_m_prime_ = std::make_unique<TransformationCompDFT>(m_prime, p, reversed_, DFT_p);

		twiddleMatrix_ = std::make_unique<entry_type_vec>(phi_m);

		// m-th rot of unity
		entry_type tempRootOfUnity = RLWE_Toolkit::Math_util::computeRootOfUnity(m);

		if (reversed_){// compute conjugated twiddle matrix
			for (int i = 1; i < p; i++)
			{
				for (int j = 0; j < m_prime; j++)
				{
					(*twiddleMatrix_)[(i - 1)*m_prime + j] = conj(pow(tempRootOfUnity, i*j));
				}
			}
		}
		else{// compute twiddle matrix
			for (int i = 1; i < p; i++)
			{
				for (int j = 0; j < m_prime; j++)
				{
					(*twiddleMatrix_)[(i - 1)*m_prime + j] = pow(tempRootOfUnity, i*j);
				}
			}
		}
	}
}

TransformationCompCRT::~TransformationCompCRT(){}

// Apply this matrix to @param vec and store result in x.
void TransformationCompCRT::applyToVector(entry_type_vec& x, entry_type_vec const& vec) const
{
	int phi_p = CRT_p_->getDim();
	int m_prime = dim_ / phi_p;
	if (m_prime <= 1){ // prime case (m = p): transformation is not decomposed. twiddleMatrix_ and DFT_p_ are null pointers.
		CRT_p_->applyToVector(x, vec);
	}
	else{ // apply matrices according to the sparse decomposition
		if (reversed_){
			applySingleKroneckerDecomposedMatrix<entry_type>(x, DFT_m_prime_.get(), vec, phi_p, 1);
			std::transform(twiddleMatrix_->cbegin(), twiddleMatrix_->cend(), x.cbegin(), x.begin(), std::multiplies<entry_type>());
			applySingleKroneckerDecomposedMatrix<entry_type>(x, CRT_p_.get(), x, 1, m_prime);
		}
		else{
			applySingleKroneckerDecomposedMatrix<entry_type>(x, CRT_p_.get(), vec, 1, m_prime);
			std::transform(twiddleMatrix_->cbegin(), twiddleMatrix_->cend(), x.cbegin(), x.begin(), std::multiplies<entry_type>());
			applySingleKroneckerDecomposedMatrix<entry_type>(x, DFT_m_prime_.get(), x, phi_p, 1);
		}
	}	
}

