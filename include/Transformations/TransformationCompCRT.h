#pragma once

#ifndef TRANS_C_CRT_H

#define TRANS_C_CRT_H

#include <tuple>

#include "Transformations/AbstractVectorTransformation.h"
#include "Transformations/MatrixCompFFT.h"

namespace RLWE_Toolkit {
	namespace Transformations {
		class TransformationCompCRT;
	}
}

/*
This class represents the Chinese remainder transformations (CRT_m, CRT_m^-1, CRT_m^*, (CRT_m^*)^-1) for prime 
powers m = p^k. We define m' = m/p. The CRT-matrix is represented by three smaller matrices: DFT_m', CRT_p and 
T_m the twiddle matrix. These matrices give a sparse decomposition which leads to a more efficient application 
of CRT.

This class is for the case where the entries are in CC, so it is used to compute the embedding sigma: K -> H.

This class implements the inherited virtual function "applyToVector". The member variables DFT_ and CRT_ are
of the same base type as this class. This->applyToVector uses applySingleKroneckerDecomposedMatrix() and therefore calls the 
applyToVector() function of DFT_ and CRT_.

@author Christoph Mayer
@version 1.0
*/
class RLWE_Toolkit::Transformations::TransformationCompCRT : public AbstractVectorTransformation<complex_type>
{
public:
	typedef AbstractVectorTransformation<complex_type> base_type;

	// indicator for the complex vector transformations
	enum class TransformationType
	{
		COMPLEX_CRT_M,
		COMPLEX_CRT_M_INVERSE,
		COMPLEX_CRT_M_STAR,
		COMPLEX_CRT_M_STAR_INVERSE,
	};

	// Constructs the by @param transformation indicated CRT matrix for prime powers m = p^k.
	// The dimension is phi(m) = (p-1)*m', where m' = m/p.
	TransformationCompCRT(int m, int p, TransformationType transformation);

	// Default destructor
	~TransformationCompCRT();


protected:

	// Apply this matrix to a given vector @param vec and store the result in x.
	virtual void applyToVector(entry_type_vec& x, entry_type_vec const& vec) const;

private:

	// Indicates if the sparse decomposition has their factors in reversed order.
	// The value is directly linked to the indicator given on construction
	bool reversed_; 

	std::unique_ptr<base_type const> CRT_p_; // actual data type: MatrixCompFFT
	std::unique_ptr<base_type const> DFT_m_prime_; // actual data type: TransformationCompDFT
	std::unique_ptr<entry_type_vec> twiddleMatrix_;	
	// diagonal "twiddle" matrix represented by a vector
};

#endif // !TRANS_C_CRT_H
