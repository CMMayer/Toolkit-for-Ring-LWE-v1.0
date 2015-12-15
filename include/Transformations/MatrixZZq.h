#pragma once

#ifndef MATRIX_ZZQ_H

#define MATRIX_ZZQ_H

#include <vector>

#include "NTL/mat_ZZ_p.h"

#include "Transformations/AbstractVectorTransformation.h"

/*
For a prime power m = p^k and m' = m/p, the matrices DFT_m',q and CRT_p,q over ZZ_q are realized by this class.

Implements the virtual function "applyToVector" and uses matrix-vector multiplication modulo q.

The matrix is realized by the NTL.

@author Christoph Mayer
@version 1.0
*/
class MatrixZZq : public AbstractVectorTransformation<integral_matrix_entry_type>
{
public:
	typedef NTL::mat_ZZ_p matrix_type;

	// Constructs an integral square matrix of size entries.NumCols() and entries from @param entries.
	MatrixZZq(matrix_type const& entries);

	// Copy constructor
	MatrixZZq(MatrixZZq const& matrix);

	// Default destructor
	~MatrixZZq();

	// Inverts this matrix modulo q. The modulus q is implicitly given by the entries of the matrix.
	// Overrides the matrix so that this class represents the inverse afterwards.
	void invert();

protected:

	// Apply this matrix to a given vector @param vec and store the result in x.
	virtual void applyToVector(entry_type_vec& x, entry_type_vec const& vec) const;

private:

	std::unique_ptr<matrix_type> matrix_;

};

#endif // !MATRIX_ZZQ_H
