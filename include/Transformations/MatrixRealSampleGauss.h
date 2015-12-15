#pragma once

#ifndef MATRIX_REAL_SAMPLE_GAUSS_H

#define MATRIX_REAL_SAMPLE_GAUSS_H

#include "boost/numeric/ublas/matrix.hpp"

#include "Transformations/AbstractVectorTransformation.h"

/*
Represents the prime-indexed real transformations C*B' used for the sampling of Gaussians in R^dual.

Implements the virtual function "applyToVector" and uses standard matrix-vector multiplication. 

@author Christoph Mayer
@version 1.0
*/
class MatrixRealSampleGauss : public AbstractVectorTransformation<real_type>
{
public:

	typedef std::vector<entry_type> entry_type_matrix;

	// Constructs C*B' of dimension @param dim
	MatrixRealSampleGauss(int dim);

	// Copy constructor
	MatrixRealSampleGauss(MatrixRealSampleGauss const& mr);

	// Default destructor
	~MatrixRealSampleGauss();

protected:

	// Apply this matrix to a given vector @param vec and store the result in x.
	virtual void applyToVector(entry_type_vec& x, entry_type_vec const& vec) const;

private:
	std::unique_ptr<entry_type_matrix> matrix_;

};

#endif // !MATRIX_REAL_SAMPLE_GAUSS_H