#pragma once

#ifndef MATRIX_REAL_GS_H

#define MATRIX_REAL_GS_H

#include "boost/numeric/ublas/matrix.hpp"

#include "Transformations/AbstractVectorTransformation.h"

namespace RLWE_Toolkit {
	namespace Transformations {
		class MatrixRealGS;
	}
}

/*
Represents the prime-indexed real transformations D_p and U_p in the Gram-Schmidt orthogonalization
of CRT_p.

Implements the virtual function "applyToVector" and uses specialized linear-time algorithms for this
multiplication.

@author Christoph Mayer
@version 1.0
*/
class RLWE_Toolkit::Transformations::MatrixRealGS : public AbstractVectorTransformation<real_type>
{
public:

	// Constructs D_p for @param D_or_U = true and U_p for @param D_or_U = false. Both matrices are of 
	// dimension @param dim. @param m_prime is needed in the application algorithm.
	MatrixRealGS(int dim, int m_prime, bool D_or_U);

	// Copy constructor
	MatrixRealGS(MatrixRealGS const& mr);

	// Default destructor
	~MatrixRealGS();

protected:

	// Apply this matrix to a given vector @param vec and store the result in x.
	virtual void applyToVector(entry_type_vec& x, entry_type_vec const& vec) const;

private:
	int m_prime_; // m' = m/p is needed in the application algorithm
	bool GS_D_or_U_; // true for D_p, false for U_p

	// Subroutines for applyToVector. GS_D_or_U_ determines which one is used.

	// Subroutine for multiplication with g in the decoding basis
	void gsU(entry_type_vec& x, entry_type_vec const& vec) const;

	// Subroutine for multiplication with g^-1 in the decoding basis
	void gsD(entry_type_vec& x, entry_type_vec const& vec) const;

};

#endif // ! MATRIX_REAL_GS_H