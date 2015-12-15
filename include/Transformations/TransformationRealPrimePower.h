#pragma once

#ifndef REAL_TRANSFORMATION_D_H

#define REAL_TRANSFORMATION_D_H


#include <vector>

#include "Transformations/AbstractVectorTransformation.h"

// forward declaration
class MatrixRealGS;
class MatrixRealSampleGauss;

typedef double real_type;

/*
Represents the prime power-indexed real transformations D_m, U_m and C*B' for a prime power m = p^k.
The general form of this transformation is A_m = A_p X I_m' (X = Kronecker product, I_m' = identity 
matrix of size m'). The dimension is phi(m) = (p-1)*m'.

This class implements the inherited virtual function "applyToVector". The member variables matrix_p_
is of the same base type as this class. This->applyToVector uses applySingleKroneckerDecomposedMatrix() 
and therefore calls the applyToVector() function of matrix_p_.

@author Christoph Mayer
@version 1.0
*/
class TransformationRealPrimePower : public AbstractVectorTransformation<real_type>
{
public:
	typedef AbstractVectorTransformation<real_type> base_type;

	// indicator for the real vector transformations
	enum class TransformationType
	{
		REAL_GS_D,				// Matrix D from Gram-Schmidt orthogonalization
		REAL_GS_U,				// Matrix U from Gram-Schmidt orthogonalization
		REAL_CONVERT_GAUSSIANS  // Matrix C*B'
	};

	// Constructs the by @param transformation indicated real transformation for prime powers m = p^k.
	// The dimension is phi(m) = (p-1)*m', where m' = m/p.
	TransformationRealPrimePower(int m, int p, TransformationType transformation);

	// Default destructor
	~TransformationRealPrimePower();

protected:

	// Apply this matrix to a given vector @param vec and store the result in x.
	virtual void applyToVector(entry_type_vec& x, entry_type_vec const& vec) const;

private:
	// Stores the prime indexed A_p from the decomposition A_m = A_p X I_m'
	std::unique_ptr<base_type const> matrix_p_; // actual data type: MatrixRealGS or MatrixRealSampleGauss

};



#endif // ! REAL_TRANSFORMATION_D_H