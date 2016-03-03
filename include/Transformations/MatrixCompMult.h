#pragma once

#ifndef MATRIX_COMP_H

#define MATRIX_COMP_H

#include "boost/numeric/ublas/matrix.hpp"

#include "Transformations/AbstractVectorTransformation.h"

namespace RLWE_Toolkit {
	namespace Transformations {
		class MatrixCompMult;
	}
}

/*
Represents the prime-indexed complex transformations CRT_p^-1 and (CRT_p^*)^-1.

Implements the virtual function "applyToVector" and uses standard matrix-vector 
multiplication. 

The matrix is realized by the boost library.

@author Christoph Mayer
@version 1.0
*/
class RLWE_Toolkit::Transformations::MatrixCompMult : public AbstractVectorTransformation<complex_type>
{
public:
	typedef boost::numeric::ublas::matrix<entry_type> matrix_type;

	// Constructs (CRT_p^*)^-1 for @param adjoint = true and CRT_p^-1 for 
	// @param adjoint = false. Both matrices are of dimension p-1.
	MatrixCompMult(int p, bool adjoint);

	// Copy constructor
	MatrixCompMult(MatrixCompMult const& matrix);

	// Default destructor
	~MatrixCompMult();

protected:

	// Apply this matrix to a given vector @param vec and store the result in x.
	virtual void applyToVector(entry_type_vec& x, entry_type_vec const& vec) const;

private:

	std::unique_ptr<matrix_type> matrix_;

	// Invert this.
	void invert();
};

#endif // !MATRIX_COMP_H