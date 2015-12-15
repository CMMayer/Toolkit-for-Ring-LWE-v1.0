#include "Transformations/MatrixZZq.h"

typedef MatrixZZq::entry_type entry_type;
typedef MatrixZZq::matrix_type matrix_type;

MatrixZZq::MatrixZZq(matrix_type const& entries)
	:
	AbstractVectorTransformation(entries.NumCols()),
	matrix_(std::make_unique<matrix_type>(entries))
{}

// Copy constructor
MatrixZZq::MatrixZZq(MatrixZZq const& entries) 
	:
	AbstractVectorTransformation(entries.getDim()),
	matrix_(std::make_unique<matrix_type>(*entries.matrix_))
{}

MatrixZZq::~MatrixZZq(){}

void MatrixZZq::invert()
{
	NTL::inv(*matrix_, *matrix_);
}

// Apply this matrix to @param vec and store result in x.
void MatrixZZq::applyToVector(entry_type_vec& x, entry_type_vec const& vec) const
{
	typedef NTL::vec_ZZ_p vector;
	vector NTL_vec = vector(NTL::INIT_SIZE, dim_);
	// make NTL_vec copy of input
	for (int i = 0; i < dim_; i++)
	{
		NTL_vec[i] = vec[i];
	}
	// multiply via NTL
	NTL_vec = *matrix_ * NTL_vec;
	for (int i = 0; i < dim_; i++)
	{
		NTL::ZZ const& r = NTL::rep(NTL_vec[i]); // rep of a coset x + pZZ of type NTL::ZZ_p is the element x of type NTL::ZZ
		x[i] = NTL::sign(r) * r.rep.rep[1]; // mult with sign to get always a positive rep. the type of r.rep.rep[1] is long.
	}
}

