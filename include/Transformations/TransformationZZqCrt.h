#pragma once

#ifndef TRANS_CRT_H

#define TRANS_CRT_H


#include <tuple>

#include "NTL/ZZ_p.h"

#include "Transformations/AbstractVectorTransformation.h"

class MatrixZZq; // forward declaration

/*
This class represents the Chinese remainder transformation (CRT) for prime powers m = p^k. We define m' = m/p. 
The CRT-matrix is represented by three smaller matrices: DFT_m', CRT_p and T_m the twiddle matrix. These matrices
give a sparse decomposition which leads to a more efficient application of CRT.

This class is for the case where the entries are in ZZ_q, q being the ring-LWE modulus. Thus it is used to switch between
the bases POWERFUL and CRT.

This class implements the inherited virtual function "applyToVector". The member variables DFT_ and CRT_ are
of the same base type as this class. This->applyToVector uses applySingleKroneckerDecomposedMatrix() and therefore calls the 
applyToVector() function of DFT_ and CRT_.

@author Christoph Mayer
@version 1.0
*/
class TransformationZZqCRT : public AbstractVectorTransformation<integral_matrix_entry_type>
{
public:
	typedef AbstractVectorTransformation<integral_matrix_entry_type> base_type;

	// Constructs the CRT_m,q matrix for prime powers m = p^k. Let m' = m/p. 
	// CRT_m is represented by three smaller matrices: DFT_m', CRT_p and the twiddle matrix.
	// @param reversed must be true, if CRT_m,q^-1 is represented, and false otherwise.
	TransformationZZqCRT(MatrixZZq const& DFT, MatrixZZq const& CRT,
		std::vector<NTL::ZZ_p> const& twiddleMatrix, bool const reversed);

	// Default destructor
	~TransformationZZqCRT();

protected:	

	// Apply this matrix to a given vector @param vec and store the result in x.
	virtual void applyToVector(entry_type_vec& x, entry_type_vec const& vec) const;

private:

	/*
	The sparse decomposition of the CRT matrix we use is only correct except for a permutation.
	See companion work [May15] for the definition of L_m^d. The value d has to be a divisor
	of m, for the definition to be well defined. m is implicitly given by a.size().
	*/
	void permutationL_M_D(entry_type_vec& x, entry_type_vec const& a, int d) const;

	// Indicates if the sparse decomposition has their factors in reversed order.
	bool reversed_;
	std::unique_ptr<base_type const> DFT_; // actual data type: MatrixZZq
	std::unique_ptr<base_type const> CRT_; // actual data type: MatrixZZq
	std::unique_ptr<std::vector<NTL::ZZ_p> const> twiddleMatrix_;	
	// diagonal matrix represented by a vector
};

#endif // !TRANS_CRT_H
