#pragma once

#ifndef TRANS_ALG_H

#define TRANS_ALG_H

#include <array>
#include <list>
#include <iostream>
#include <algorithm>
#include <functional>

#include "Transformations/AbstractVectorTransformation.h"

template<typename T> using Trans = AbstractVectorTransformation<T>;
template<typename T> using avt_uptr = std::unique_ptr < Trans<T> >;

/*
Let A = @param matrix. Then we have x = (I_dim1 X A X I_dim2) * vec. (X = Kronecker product, I_dim = identity matrix of size dim)
Computes the right sub vectors x' for multiplication with A.

@param x: Vector where the result is stored.
@param matrix: The transformation matrix.
@param vec: Vector which the transformation is applied to.
@param dim1: Dimension of the left-hand identity matrix.
@param dim2: Dimension of the right-hand identity matrix.
*/
template<typename T> void applySingleKroneckerDecomposedMatrix(std::vector<T>& x, Trans<T> const* matrix,
	std::vector<T> const& vec, int const dim1, int const dim2)
{
	int n = vec.size();
	int const m = matrix->getDim();
	std::vector<T> apply_vec(m, 0); //vector for applying the matrix, apply_vec = matrix * apply_vec

	for (int u = 0; u < n; u += (n / dim1))
	{
		for (int v = 0; v < dim2; v++)
		{
			for (int w = 0; w < m; w++)
			{ // compute appropriate sub vector 
				apply_vec[w] = vec[w*dim2 + v + u]; 
			}
			try{ // and apply the matrix to this sub vector
				matrix->applyToVector(apply_vec, apply_vec);
			}
			catch (...){ // the decoding procedure might throw an exception
				throw;
			}
			for (int w = 0; w < m; w++)
			{ // update the resulting x vector
				x[w*dim2 + v + u] = apply_vec[w];
			}
		}
	}
}

/*
This procedure provides an efficient computation for the application of Kronecker decomposed matrices. 
The given list represent a Kronecker decomposed matrix A = A_1 X ... X A_s. The algorithms produces
x = A * vec.
Calls applySingleKroneckerDecomposedMatrix().

@param x: Vector where the result is stored.
@param matrices: List of transformation matrices representing the Kronecker decomposition.
@param vec: Vector which the transformation is applied to.
*/
template<typename T> void applyKroneckerDecomposition(std::vector<T>& x, std::list<avt_uptr<T>> const& matrices,
	std::vector<T> const& vec)
{
	int const n = vec.size();
	int const m = matrices.size();	//number of matrices in the Kronecker decomposition.
	std::vector<int> matrix_sizes(m);

	if (!(x == vec)){	// make x a copy of vec if there not equal. We use x for temporary computations and the output.
		x = vec;
	}
	
	// store the smaller matrix dimensions in a vector.
	std::transform(matrices.cbegin(), matrices.cend(), matrix_sizes.begin(), [](avt_uptr<T> const& mat){ return mat->getDim(); });
	
	int i = 0; // separate counter for the "for each" loop.
	
	for each (auto const& matrix in matrices)
	{
		int dim1 = 1, dim2 = 1; // initiate as the multiplicative neutral element.

		// compute dimension of the left-hand identity matrix
		for (int k = 0; k < i; k++)
		{
			dim1 *= matrix_sizes[k];
		}

		// compute dimension of the right-hand identity matrix
		for (int k = i + 1; k < m; k++)
		{
			dim2 *= matrix_sizes[k];
		}
		try{
			applySingleKroneckerDecomposedMatrix<T>(x, matrix.get(), x, dim1, dim2);
		}
		catch (...){
			throw;
		}
		i++;
	}
}

#endif // !TRANS_ALG_H