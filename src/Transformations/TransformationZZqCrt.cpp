#include "Transformations/TransformationZZqCRT.h"
#include "Transformations/Transformation algorithms.h"
#include "Transformations/MatrixZZq.h"

typedef TransformationZZqCRT::entry_type entry_type;

TransformationZZqCRT::TransformationZZqCRT(MatrixZZq const& DFT, MatrixZZq const& CRT, std::vector<NTL::ZZ_p> const& twiddleMatrix,
	bool const reversed)
	:
	AbstractVectorTransformation(twiddleMatrix.size()),
	DFT_(std::make_unique<MatrixZZq>(DFT)),
	CRT_(std::make_unique<MatrixZZq>(CRT)),
	twiddleMatrix_(std::make_unique<std::vector<NTL::ZZ_p>>(twiddleMatrix)),
	reversed_(reversed)
{}

TransformationZZqCRT::~TransformationZZqCRT(){}

// Apply this matrix to @param vec and store result in x.
void TransformationZZqCRT::applyToVector(entry_type_vec& x, entry_type_vec const& vec)const
{
	int phi_p = CRT_->getDim();
	int m_prime = DFT_->getDim();

	// local function for multiplication of integers with NTL::ZZ_p objects
	std::function<int(NTL::ZZ_p, entry_type)> multZZp =
		[](NTL::ZZ_p x, entry_type y)
	{
		NTL::ZZ_p y2 = NTL::ZZ_p(y);
		y2 *= x;
		return NTL::rep(y2).rep.rep[1];
	};

	// apply matrices according to the sparse decomposition
	if (reversed_){
		//In the discrete case CRT_m,q, the permutation is important (in contrast to the complex case CRT_m)
		permutationL_M_D(x, vec, phi_p); 

		applySingleKroneckerDecomposedMatrix<entry_type>(x, DFT_.get(), x, phi_p, 1);
		std::transform(twiddleMatrix_->cbegin(), twiddleMatrix_->cend(), x.cbegin(), x.begin(), multZZp);
		applySingleKroneckerDecomposedMatrix<entry_type>(x, CRT_.get(), x, 1, m_prime);
	} 
	else{
		applySingleKroneckerDecomposedMatrix<entry_type>(x, CRT_.get(), vec, 1, m_prime);
		std::transform(twiddleMatrix_->cbegin(), twiddleMatrix_->cend(), x.cbegin(), x.begin(), multZZp);
		applySingleKroneckerDecomposedMatrix<entry_type>(x, DFT_.get(), x, phi_p, 1);

		//In the discrete case CRT_m,q, the permutation is important (in contrast to the complex case CRT_m)
		permutationL_M_D(x, x, m_prime);
	}
}

// The permutation L_m^d. d has to be a divisor of m = a.size(). See [May15] for a definition of L_m^d.
void TransformationZZqCRT::permutationL_M_D(entry_type_vec& x, entry_type_vec const& a, int d) const
{
	int m = a.size();
	entry_type_vec res = entry_type_vec(m);

	res[0] = a[0];
	res[m - 1] = a[m - 1];
	for (int i = 1; i < m - 1; i++)
	{
		res[i] = a[i*d % (m - 1)];
	}
	x = res;
}

