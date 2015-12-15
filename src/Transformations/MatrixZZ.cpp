#include "Transformations/MatrixZZ.h"
#include <algorithm>

typedef MatrixZZ::MatrixType MatrixType;

MatrixZZ::MatrixZZ(int dim, MatrixType matrixType)
	:
	AbstractVectorTransformation(dim),
	matrixType_(matrixType)
{}

// Copy constructor
MatrixZZ::MatrixZZ(MatrixZZ const& mzz)
	:
	AbstractVectorTransformation(mzz.dim_),
	matrixType_(mzz.getMultType())
{}

MatrixZZ::~MatrixZZ(){}

// Apply this matrix to @param vec and store result in x.
void MatrixZZ::applyToVector(entry_type_vec& x, entry_type_vec const& vec) const
{
	switch (matrixType_)
	{ // Decide which specialized application function is used.
	case MatrixType::MATRIX_L:
		applyL(x, vec);
		break;

	case MatrixType::MATRIX_L_INVERSE:
		applyLInverse(x, vec);
		break;

	case MatrixType::MATRIX_G_DECODING:
		gDec(x, vec);
		break;

	case MatrixType::MATRIX_G_INVERSE_DECODING:
		try{
			gInvDec(x, vec);
		}
		catch (std::exception& e){ // a decoding error might occur
			throw e;
		}
		break;

	case MatrixType::MATRIX_G_POWERFUL:
		x = vec;
		std::reverse(x.begin(), x.end());// Corresponds to multiplication with reverse identity matrix J.
		gPow(x, x);
		std::reverse(x.begin(), x.end());// Corresponds to multiplication with reverse identity matrix J.
		break;

	case MatrixType::MATRIX_G_INVERSE_POWERFUL:
		try{
			x = vec;
			std::reverse(x.begin(), x.end());// Corresponds to multiplication with reverse identity matrix J.
			gInvPow(x, x);
			std::reverse(x.begin(), x.end());// Corresponds to multiplication with reverse identity matrix J.
		}
		catch (std::exception& e){ // a decoding error might occur
			throw e;
		}
		break;
	}
}

// Specialized application function for G_p^dec. See definition of G_p^dec in [May15]
void MatrixZZ::gDec(entry_type_vec& x, entry_type_vec const& vec) const
{
	int n = dim_;
	entry_type res = 0;
	for (int i = n - 1; i > 0; i--)
	{
		res += vec[i];
		x[i] = vec[i] - vec[i - 1];
	}
	x[0] = res + 2 * vec[0];
}

// Specialized application function for (G_p^dec)^-1. See definition of (G_p^dec)^-1 in [May15]
void MatrixZZ::gInvDec(entry_type_vec& x, entry_type_vec const& vec) const
{
	int n = dim_;
	int p = dim_ + 1;
	entry_type res = vec[0];
	for (int i = 1; i < n; i++)
	{
		res += (i + 1 - p)*vec[i];
	}
	x[0] = res;
	for (int i = 1; i < n; i++)
	{
		res += p*vec[i];
		x[i] = res;
	}
	try{
		std::for_each(x.begin(), x.end(), [p](entry_type& e){
			if (e % p == 0){
				e /= p;
			}
			else{
				throw RingLweException("Couldn't divide by p while decoding procedure. Decoding failure detected!\n\n");
			}
		});
	}
	catch (std::exception& e){
		throw e;
	}
}

// Specialized application function for G_p^pow. See definition of G_p^pow in [May15]
void MatrixZZ::gPow(entry_type_vec& x, entry_type_vec const& vec) const
{
	int n = dim_;
	entry_type_vec temp = entry_type_vec(n);
	temp[0] = 2 * vec[0] - vec[1];
	for (int i = 1; i < n-1; i++)
	{
		temp[i] = vec[0] + vec[i] - vec[i + 1];
	}
	temp[n - 1] = vec[0] + vec[n - 1];
	x = temp;
}

// Specialized application function for (G_p^pow)^-1. See definition of (G_p^pow)^-1 in [May15]
void MatrixZZ::gInvPow(entry_type_vec& x, entry_type_vec const& vec) const
{
	int n = dim_;
	int p = dim_ + 1;
	entry_type_vec temp = entry_type_vec(n);
	entry_type res = 0;
	for (int i = 0; i < n; i++)
	{
		res += vec[i];
	}
	entry_type res2 = res;
	for (int i = 0; i < n; i++)
	{
		temp[i] = res2;
		res2 += res - p*vec[i];
	}
	x = temp;
	std::for_each(x.begin(), x.end(), [p](entry_type& e){
		if (e % p == 0){
			e /= p;
		}
		else{
			throw RingLweException("Couldn't divide by p while decoding procedure. Decoding failure detected!\n\n");
		}
	});
}

// Specialized application function for L_p. See definition of L_p in [May15]
void MatrixZZ::applyL(entry_type_vec& x, entry_type_vec const& vec) const
{
	entry_type res = 0;
	for (int i = 0; i < dim_; i++)
	{
		res += vec[i];
		x[i] = res;
	}
}

// Specialized application function for L_p^-1. See definition of L_p^-1 in [May15]
void MatrixZZ::applyLInverse(entry_type_vec& x, entry_type_vec const& vec) const
{
	entry_type temp1, temp2 = vec[0];
	x[0] = temp2;
	for (int i = 1; i < dim_; i++)
	{
		temp1 = temp2;
		temp2 = vec[i];
		x[i] = vec[i] - temp1;
	}
}

