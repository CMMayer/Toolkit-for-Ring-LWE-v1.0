#pragma once

#ifndef TRANS_C_DFT_H

#define TRANS_C_DFT_H

#include <tuple>

#include "Transformations/AbstractVectorTransformation.h"
#include "Transformations/MatrixCompFFT.h"

/*
This class represents the discrete Fourier transformations (DFT_m, DFT_m^-1, DFT_m^*, (DFT_m^*)^-1) for prime
powers m = p^k. We define m' = m/p. The DFT-matrix is represented by three smaller matrices: DFT_m', DFT_p and 
T_m the twiddle matrix. These matrices give a sparse decomposition which leads to a more efficient application 
of DFT.

This class is for the case where the entries are in CC, so it is used to compute the embedding sigma: K -> H.

This class implements the inherited virtual function "applyToVector". The member variables DFT_m_prime and DFT_p 
are of the same base type as this class. This->applyToVector uses applySingleKroneckerDecomposedMatrix() and therefore calls the
applyToVector() function of DFT_m_prime and DFT_p.

@author Christoph Mayer
@version 1.0
*/
class TransformationCompDFT : public AbstractVectorTransformation<complex_type>
{
public:
	typedef AbstractVectorTransformation<complex_type> base_type;

	// Constructs the DFT matrix for prime powers m = p^k. Calls itself recursively until m = p. 
	// @param reversed must be true, if DFT_m^-1 or DFT_m^* is represented, and false otherwise.
	// @param DFT_p is the respective prime indexed DFT in the sparse decomposition.
	TransformationCompDFT(int m, int p, bool reversed, std::shared_ptr<MatrixCompFFT const> DFT_p);

	// Default destructor
	~TransformationCompDFT();


protected:

	// Apply this matrix to a given vector @param vec and store the result in x.
	virtual void applyToVector(entry_type_vec& x, entry_type_vec const& vec) const;

private:

	// Indicates if the sparse decomposition has their factors in reversed order.
	bool reversed_;

	std::shared_ptr<base_type const> DFT_p_; // actual data type: MatrixCompFFT
	std::unique_ptr<base_type const> DFT_m_prime_; // actual data type: TransformationCompDFT
	std::unique_ptr<entry_type_vec> twiddleMatrix_;	
	// diagonal "twiddle" matrix represented by a vector
};

#endif // !TRANS_C_DFT_H
