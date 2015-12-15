#pragma once

#ifndef MATRIX_MG_H

#define MATRIX_MG_H

#include <vector>

#include "Transformations/AbstractVectorTransformation.h"

/*
Represents the prime-indexed integral vector transformations G_p^dec, G_p^pow and their inverses
for multiplication with the special elements g and g^-1 in the powerful and decoding bases as well 
as the integral matrices L_p and L_p^-1 for conversion between t^-1 p and d.

Implements the virtual function "applyToVector" and uses specialized linear-time algorithms for this
multiplication.

@author Christoph Mayer
@version 1.0
*/
class MatrixZZ : public AbstractVectorTransformation<integral_matrix_entry_type>
{
public:

	// Indicator for the different represented integer transformations
	enum class MatrixType
	{
		MATRIX_L,				// the matrix L
		MATRIX_L_INVERSE,		// the matrix L^-1
		MATRIX_G_DECODING,				// represents multiplication with g in decoding
		MATRIX_G_INVERSE_DECODING,		// represents multiplication with g^-1 in decoding
		MATRIX_G_POWERFUL,				// represents multiplication with g in powerful
		MATRIX_G_INVERSE_POWERFUL		// represents multiplication with g^-1 in powerful
	};

	// Constructs the integer vector transformation indicated by @param matrixType.
	// The dimension of the transformation is @param dim.
	MatrixZZ(int dim, MatrixType matrixType);

	// Copy constructor
	MatrixZZ(MatrixZZ const& mzz);

	// Default destructor
	~MatrixZZ();

	inline MatrixType getMultType() const { return matrixType_; };

protected:

	// Apply this matrix to a given vector @param vec and store the result in x.
	virtual void applyToVector(entry_type_vec& x, entry_type_vec const& vec) const;

private:

	MatrixType matrixType_;

	// Subroutines for applyToVector. matrixType_ determines which one is used.

	// Subroutine for multiplication with g in the decoding basis
	void gDec(entry_type_vec& x, entry_type_vec const& vec) const;

	// Subroutine for multiplication with g^-1 in the decoding basis
	void gInvDec(entry_type_vec& x, entry_type_vec const& vec) const;

	// Subroutine for multiplication with g in the powerful basis
	void gPow(entry_type_vec& x, entry_type_vec const& vec) const;

	// Subroutine for multiplication with g^-1 in the powerful basis
	void gInvPow(entry_type_vec& x, entry_type_vec const& vec) const;

	// Subroutine for application of L
	void applyL(entry_type_vec& x, entry_type_vec const& vec) const;

	// Subroutine for application of L^-1
	void applyLInverse(entry_type_vec& x, entry_type_vec const& vec) const;
};

#endif // ! MATRIX_MG_H