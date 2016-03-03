#pragma once

#ifndef MATRIX_COMP_FFT_H

#define MATRIX_COMP_FFT_H

#include "Transformations/AbstractVectorTransformation.h"

namespace RLWE_Toolkit {
	namespace Transformations {
		class MatrixCompFFT;
	}
}

/*
Represents the prime-indexed complex transformations, e.g., DFT_p or CRT_p. 
(An exception is made by CRT_p^-1 and (CRT_p^*)^-1

Implements the virtual function "applyToVector" and uses FFT algorithms for this 
multiplication.

@author Christoph Mayer
@version 1.0
*/
class RLWE_Toolkit::Transformations::MatrixCompFFT : public AbstractVectorTransformation<complex_type>
{
public:
	typedef std::vector<entry_type> entry_type_vec;

	// Indicator for the different complex transformations represented by this class
	enum class MatrixType{
		DFT_P,	
		DFT_P_INV,		// DFT_p^-1
		DFT_P_STAR,		// DFT_p^*
		DFT_P_STAR_INV, // (DFT_p^*)^-1
		CRT_P,
		CRT_P_STAR		// CRT_p^*
	};

	// Constructs the by @param trans indicated transformation of dimension p or p-1.
	MatrixCompFFT(MatrixType matrixType, int p);

	~MatrixCompFFT();

protected:

	// Apply this matrix to a given vector @param vec and store the result in x.
	virtual void applyToVector(entry_type_vec& x, entry_type_vec const& vec) const;

private:

	MatrixType matrixType_;
	int generator; // A generator of the group ZZ_P^*
	std::unique_ptr<entry_type_vec> precomp_DFT_omega_p_;	// Pre-computed DFT of powers of the p-th root of unity

};

#endif // !MATRIX_COMP_FFT_H
