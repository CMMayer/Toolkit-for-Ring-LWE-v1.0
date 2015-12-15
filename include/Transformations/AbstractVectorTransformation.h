#pragma once

#ifndef ABS_VEC_TRANS_H

#define ABS_VEC_TRANS_H

#include <vector>
#include <list>
#include <memory>
#include <iostream>
#include <complex>
#include "Exceptions/RingLweException.h"

typedef int integral_matrix_entry_type;  // entry_type for all integral matrices (int or long)
typedef double real_type; // type for representation of the reals. Also for usage in std::complex. 
typedef std::complex<real_type> complex_type;

/*
The abstract base class for all square vector transformations. Typename T is the entry type for all matrix and vector entries.
Provides pure virtual function "applyToVector".

@author Christoph Mayer
@version 1.0
*/
template<typename T> class AbstractVectorTransformation
{
public:
	typedef T entry_type;
	typedef std::vector<entry_type> entry_type_vec;

	// Getter
	inline int getDim() const { return dim_; };

	// Purely virtual function to apply the matrix to a given vector @param vec and store the result in x.
	virtual void applyToVector(entry_type_vec& x, entry_type_vec const& vec) const = 0;

	~AbstractVectorTransformation(){};

protected:

	AbstractVectorTransformation(int dim):dim_(dim){};

	int dim_;	// Dimension of the transformation matrix.
};

#endif // !ABS_VEC_TRANS_H