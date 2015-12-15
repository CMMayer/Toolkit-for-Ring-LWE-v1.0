#pragma once

#ifndef TRANS_LM_H

#define TRANS_LM_H


#include <vector>

#include "Transformations/AbstractVectorTransformation.h"

class MatrixZZ; // forward declaration

/*
Represents the prime power-indexed integer transformations L_m, G_m^dec, G_m^pow and their inverses 
for a prime power m = p^k. The general form of this transformation is A_m = A_p X I_m' (X = Kronecker 
product, I_m' = identity matrix of size m'). The dimension is phi(m) = (p-1)*m'.

This class implements the inherited virtual function "applyToVector". The member variables matrix_p_
is of the same base type as this class. This->applyToVector uses applySingleKroneckerDecomposedMatrix()
and therefore calls the applyToVector() function of matrix_p_.

@author Christoph Mayer
@version 1.0
*/
class TransformationZZ_PrimePower : public AbstractVectorTransformation<integral_matrix_entry_type>
{
public:
	typedef AbstractVectorTransformation<integral_matrix_entry_type> base_type;

	// indicator for the integer vector transformations
	enum class TransformationType
	{
		INT_L,						// the matrix L
		INT_L_INVERSE,				// the matrix L^-1
		INT_G_DECODING,				// represents multiplication with g in decoding
		INT_G_INVERSE_DECODING,		// represents multiplication with g^-1 in decoding
		INT_G_POWERFUL,				// represents multiplication with g in powerful
		INT_G_INVERSE_POWERFUL		// represents multiplication with g^-1 in powerful
	};

	// Constructs the by @param transformation indicated integer transformation for prime powers m = p^k.
	// The dimension is phi(m) = (p-1)*m', where m' = m/p.
	TransformationZZ_PrimePower(int m, int p, TransformationType transformation);

	// Default destructor
	~TransformationZZ_PrimePower();

protected:

	// Apply this matrix to a given vector @param vec and store the result in x.
	virtual void applyToVector(entry_type_vec& x, entry_type_vec const& vec) const;

private:
	// Stores the prime indexed A_p from the decomposition A_m = A_p X I_m'
	std::unique_ptr<base_type const> matrix_p_; // actual data type: MatrixZZ

};

#endif // !TRANS_LM_H