#define _USE_MATH_DEFINES

#include <math.h>

#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/vector_proxy.hpp"
#include "boost/numeric/ublas/triangular.hpp"
#include "boost/numeric/ublas/lu.hpp"
#include "boost/numeric/ublas/io.hpp"

#include "Transformations/MatrixCompMult.h"
#include "Math util.h"

using namespace boost::numeric::ublas;

using RLWE_Toolkit::Transformations::MatrixCompMult;

typedef MatrixCompMult::entry_type entry_type;
typedef MatrixCompMult::matrix_type matrix_type;

MatrixCompMult::MatrixCompMult(int p, bool adjoint) :
	AbstractVectorTransformation(p - 1)
{
	matrix_ = std::make_unique<matrix_type>(dim_, dim_);

	// p-th root of unity
	entry_type tempRootOfUnity = RLWE_Toolkit::Math_util::computeRootOfUnity(p);

	// compute CRT for prime p with dim(CRT_p) = phi(p) = p-1.
	for (pos_int i = 1; i < p; i++)
	{
		for (pos_int j = 0; j < p - 1; j++)
		{
			(*matrix_)(i - 1, j) = pow(tempRootOfUnity, i*j);
		}
	}
	invert(); // compute CRT_m^-1
	if (adjoint){ // compute (CRT_m^*)^-1
		*matrix_ = herm(*matrix_);
	}
}

MatrixCompMult::MatrixCompMult(MatrixCompMult const& entries) :
	AbstractVectorTransformation(entries.getDim()),
	matrix_(std::make_unique<matrix_type>(*entries.matrix_))
{}

MatrixCompMult::~MatrixCompMult(){}

// Apply this matrix to @param vec and store result in x.
void MatrixCompMult::applyToVector(entry_type_vec& x, entry_type_vec const& vec)const
{
	typedef boost::numeric::ublas::vector<entry_type, entry_type_vec> vector;
	vector boost_vec = vector(vec.size());
	// make boost_vec copy of input
	std::copy(vec.cbegin(), vec.cend(), boost_vec.begin());
	// multiply via boost library
	boost_vec = prod(*matrix_, boost_vec);
	// extract result
	x = boost_vec.data();
}

// invert via LU-factorization
void MatrixCompMult::invert()
{
	using namespace boost::numeric::ublas;
	typedef permutation_matrix<std::size_t> pmatrix;
	// create a working copy of the input
	matrix<entry_type> A(*matrix_);
	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);

	// create identity matrix of "inverse"
	matrix_->assign(identity_matrix<entry_type>(A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, *matrix_);
}

