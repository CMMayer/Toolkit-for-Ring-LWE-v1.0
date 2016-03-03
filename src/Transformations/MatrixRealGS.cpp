#include "Transformations/MatrixRealGS.h"
#include <algorithm>

using RLWE_Toolkit::Transformations::MatrixRealGS;

MatrixRealGS::MatrixRealGS(int dim, int m_prime, bool D_or_U) :
	AbstractVectorTransformation(dim),
	m_prime_(m_prime),
	GS_D_or_U_(D_or_U)
{}

// Copy constructor
MatrixRealGS::MatrixRealGS(MatrixRealGS const& mr) :
	AbstractVectorTransformation(mr.dim_),
	m_prime_(mr.m_prime_),
	GS_D_or_U_(mr.GS_D_or_U_)
{}

MatrixRealGS::~MatrixRealGS(){}

// Apply this matrix to @param vec and store result in x.
void MatrixRealGS::applyToVector(entry_type_vec& x, entry_type_vec const& vec) const
{
	if (GS_D_or_U_)
	{
		gsD(x, vec);
	}
	else{
		gsU(x, vec);
	}
}

// Specialized application function for D_p. See definition of D_p in [May15]
void MatrixRealGS::gsD(entry_type_vec& x, entry_type_vec const& vec) const
{
	int p = dim_ + 1;
	for (int i = 0; i < dim_; i++)
	{
		x[i] = sqrt(m_prime_ * ((p - 1) - i / (p - i))) * vec[i];
	}
}

// Specialized application function for U_p. See definition of U_p in [May15]
void MatrixRealGS::gsU(entry_type_vec& x, entry_type_vec const& vec) const
{
	int p = dim_ + 1;
	entry_type_vec temp = vec;
	entry_type res = temp[0];
	x[0] = res;
	for (int i = 1; i < dim_; i++)
	{
		res = res + ((-1 / (p - i) - 1)*temp[i - 1] + temp[i]);
		x[i] = res;
	}
}

