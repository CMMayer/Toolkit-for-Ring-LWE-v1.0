#define _USE_MATH_DEFINES


#include <math.h>
#include <random>
#include <iostream>
#include <numeric>

#include "NTL/mat_ZZ_p.h"
#include "NTL/vec_ZZ_p.h"
#include "NTL/ZZ.h"

#include "Math util.h"

#include "Main_Classes/RingLweCryptographyField.h"
#include "Main_Classes/RingLweCryptographyElement.h"
#include "Transformations/TransformationCompCRT.h"
#include "Transformations/TransformationZZqCRT.h"
#include "Transformations/TransformationZZ_PrimePower.h"
#include "Transformations/TransformationRealPrimePower.h"
#include "Transformations/MatrixCompMult.h"
#include "Transformations/MatrixZZ.h"
#include "Transformations/MatrixZZq.h"
#include "Transformations/MatrixRealGS.h"
#include "Transformations/MatrixRealSampleGauss.h"

typedef RLWE_Toolkit::Main_Classes::RingLweCryptographyElement RingElement;
typedef RLWE_Toolkit::Main_Classes::RingLweCryptographyField LweField;
typedef RingElement::Basis Basis;
typedef LweField::TransformationMatrices TransMat;
 

LweField::RingLweCryptographyField(pos_int m, pos_int modulus) :
m_(m),
n_(RLWE_Toolkit::Math_util::eulerTotient(m_)),
rad_m_(1),
rlwe_modulus_(modulus)
{
	init();
}

LweField::~RingLweCryptographyField(){}

// Returns the 0 with basis @param basis
RingElement LweField::getZero(RingElement::Basis basis) const
{
	coordinate_vec coords(getDimension(), 0);
	RingElement x(shared_from_this(), basis, coords);
	if (basis == Basis::BASIS_CRT || basis == Basis::BASIS_T_INVERSE_CRT){
		x.setPrimeModulusQ(rlwe_modulus_);
	}
	return x;
}

// Returns the 1 with basis @param basis
RingElement LweField::getOne(RingElement::Basis basis) const
{
	RingElement x(shared_from_this(), basis, *(*ones)[basis]);
	if (basis == Basis::BASIS_CRT){
		x.setPrimeModulusQ(rlwe_modulus_);
	}
	return x;
}

void LweField::applyIntTransformation(TransformationMatrices transformation, coordinate_vec& x, coordinate_vec const& a) const
{
	switch (transformation)
	{
	case TransMat::INTEGER_CRT_MQ:
		applyKroneckerDecomposition<coordinate_type>(x, *crtMq_, a);
		return;
	case TransMat::INTEGER_CRT_MQ_INVERSE:
		applyKroneckerDecomposition<coordinate_type>(x, *crtMqInverse_, a);
		return;
	case TransMat::INTEGER_L_M:
		applyKroneckerDecomposition<coordinate_type>(x, *Lm_, a);
		return;
	case TransMat::INTEGER_L_M_INVERSE:
		applyKroneckerDecomposition<coordinate_type>(x, *LmInverse_, a);
		return;
	case TransMat::INTEGER_MULT_G_DEC:
		applyKroneckerDecomposition<coordinate_type>(x, *mgd_, a);
		return;
	case TransMat::INTEGER_DIV_G_DEC:
		try{
			applyKroneckerDecomposition<coordinate_type>(x, *dgd_, a);
		}
		catch (...){ // decoding failure might occur
			throw;
		}		
		return;
	case TransMat::INTEGER_MULT_G_POW:
		applyKroneckerDecomposition<coordinate_type>(x, *mgp_, a);
		return;
	case TransMat::INTEGER_DIV_G_POW:
		try{
			applyKroneckerDecomposition<coordinate_type>(x, *dgp_, a);
		}
		catch (...){ // decoding failure might occur
			throw;
		}
		return;
	
	default:
		return;
	}
}

void LweField::applyRealTransformation(TransformationMatrices transformation, real_vec& x, real_vec const& a) const
{
	switch (transformation)
	{
	case TransMat::REAL_SAMPLE_GAUSS_D:
		applyKroneckerDecomposition<real_type>(x, *sample_gauss_D_, a);
		return;
	case TransMat::REAL_GS_DECOMP_U:
		applyKroneckerDecomposition<real_type>(x, *gs_decomp_U_, a);
		return;
	case TransMat::REAL_GS_DECOMP_D:
		applyKroneckerDecomposition<real_type>(x, *gs_decomp_D_, a);
		return;
	default:
		return;
	}
}

void LweField::applyCompTransformation(TransformationMatrices transformation, comp_vec& x, comp_vec const& a) const
{
	switch (transformation)
	{
	case TransMat::COMPLEX_CRT_M:
		applyKroneckerDecomposition<complex_type>(x, *crtM_, a);
		return;
	case TransMat::COMPLEX_CRT_M_STAR_INVERSE:
		applyKroneckerDecomposition<complex_type>(x, *crtMStarInverse_, a);
		return;
	case TransMat::COMPLEX_CRT_M_INVERSE:
		applyKroneckerDecomposition<complex_type>(x, *crtMInverse_, a);
		return;
	case TransMat::COMPLEX_CRT_M_STAR:
		applyKroneckerDecomposition<complex_type>(x, *crtMStar_, a);
		return;

	default:
		return;
	}
}

void LweField::init()
{
	std::list<std::pair<pos_int, pos_int>> primeFactors = RLWE_Toolkit::Math_util::primeFactorization(m_);

	// Initiate all vector transformation lists
	typedef std::list<avt_uptr<coordinate_type>> uptr_type_i;
	crtMq_ = std::make_unique<uptr_type_i>();
	crtMqInverse_ = std::make_unique<uptr_type_i>();
	Lm_ = std::make_unique<uptr_type_i>();
	LmInverse_ = std::make_unique<uptr_type_i>();
	mgd_ = std::make_unique<uptr_type_i>();
	dgd_ = std::make_unique<uptr_type_i>();
	mgp_ = std::make_unique<uptr_type_i>();
	dgp_ = std::make_unique<uptr_type_i>();

	typedef std::list<avt_uptr<complex_type>> uptr_type_c;
	crtM_ = std::make_unique<uptr_type_c>();
	crtMInverse_ = std::make_unique<uptr_type_c>();
	crtMStar_ = std::make_unique<uptr_type_c>();
	crtMStarInverse_ = std::make_unique<uptr_type_c>();

	typedef std::list<avt_uptr<real_type>> uptr_type_r;
	sample_gauss_D_ = std::make_unique<uptr_type_r>();
	gs_decomp_U_ = std::make_unique<uptr_type_r>();
	gs_decomp_D_ = std::make_unique<uptr_type_r>();

	for each (std::pair<pos_int, pos_int> primePower in primeFactors){
		pos_int p = primePower.first;
		pos_int m_prime = pow(primePower.first, primePower.second - 1); // p^k-1
		pos_int m = p * m_prime;
		pos_int phi_m = (p - 1) * m_prime; // phi(p)*m' = phi(m)

		rad_m_ *= p;

		computeIntegerTransMatrices(p, m_prime, m, phi_m);
		computeComplexTransMatrices(p, m);
		computeRealTransMatrices(p, m);
	}

	computeOnes();
	computeGinCRT(primeFactors);
	computeT_inverseInD();
}

// Pre-computes all integral vector transformation matrices.
void LweField::computeIntegerTransMatrices(pos_int p, pos_int m_prime, pos_int m, pos_int phi_m)
{
	using namespace NTL;
	using namespace RLWE_Toolkit::Transformations;
	using namespace RLWE_Toolkit::Math_util;

	// init the r-LWE modulus q for the elements over ZZ_q
	ZZ q = ZZ(rlwe_modulus_);
	ZZ_p::init(q);

	//--------------------------compute CRT_p, DFT_m and twiddleMatrix----------------------------\\

	mat_ZZ_p CRT = mat_ZZ_p();
	CRT.SetDims(p - 1, p - 1);
	mat_ZZ_p DFT = mat_ZZ_p();
	DFT.SetDims(m_prime, m_prime);

	std::vector<ZZ_p> twiddleMatrix = std::vector<ZZ_p>(phi_m);

	ZZ_p integralRootOfUnity = findElementOfOrder(p, rlwe_modulus_);

	for (pos_int i = 1; i < p; i++)
	{
		for (pos_int j = 0; j < p - 1; j++)
		{
			CRT[i - 1][j] = power(integralRootOfUnity, (i*j));
		}
	}

	integralRootOfUnity = findElementOfOrder(m_prime, rlwe_modulus_);

	for (pos_int i = 0; i < m_prime; i++)
	{
		for (pos_int j = 0; j < m_prime; j++)
		{
			DFT[i][j] = power(integralRootOfUnity, (i*j));
		}
	}

	integralRootOfUnity = findElementOfOrder(m, rlwe_modulus_);

	for (pos_int i = 1; i < p; i++)
	{
		for (pos_int j = 0; j < m_prime; j++)
		{
			twiddleMatrix[(i - 1)*m_prime + j] = power(integralRootOfUnity, (i*j));
		}
	}


	//---------------------------create and push CRT_MQ------------------------------------------\\

	MatrixZZq CRT_matrix = MatrixZZq(CRT);
	MatrixZZq DFT_matrix = MatrixZZq(DFT);

	std::unique_ptr<TransformationZZqCRT::base_type> KronDecompFacCrt = std::make_unique<TransformationZZqCRT>(DFT_matrix, CRT_matrix, twiddleMatrix,
		false);

	crtMq_->push_back(std::move(KronDecompFacCrt));


	//---------------------------create and push CRT_MQ_INVERSE------------------------------------------\\

	CRT_matrix.invert();
	DFT_matrix.invert();

	std::for_each(twiddleMatrix.begin(), twiddleMatrix.end(), [](ZZ_p& x){ NTL::inv(x, x); });

	KronDecompFacCrt = std::make_unique<TransformationZZqCRT>(DFT_matrix, CRT_matrix, twiddleMatrix,
		true);

	crtMqInverse_->push_back(std::move(KronDecompFacCrt));


	// create all integral transformations with a matrix M of size phi(p) = p - 1 and actual transformation
	// size n = phi(m) (m is a prime power) as the transformation is M x I_m (Kronecker product).

	typedef TransformationZZ_PrimePower::TransformationType TransformationType;

	// one factor of the Kronecker decomposition of L_m for arbitrary m.
	std::unique_ptr<TransformationZZ_PrimePower::base_type> KronDecompFac =
		std::make_unique<TransformationZZ_PrimePower>(m, p, TransformationType::INT_L);
	Lm_->push_back(std::move(KronDecompFac));

	// one factor of the Kronecker decomposition of L_m^-1 for arbitrary m.
	KronDecompFac =
		std::make_unique<TransformationZZ_PrimePower>(m, p, TransformationType::INT_L_INVERSE);
	LmInverse_->push_back(std::move(KronDecompFac));

	// one factor of the Kronecker decomposition of G_m^d for arbitrary m.
	KronDecompFac = 
		std::make_unique<TransformationZZ_PrimePower>(m, p, TransformationType::INT_G_DECODING);
	mgd_->push_back(std::move(KronDecompFac));

	// one factor of the Kronecker decomposition of (G_m^d)^-1 for arbitrary m.
	KronDecompFac =
		std::make_unique<TransformationZZ_PrimePower>(m, p, TransformationType::INT_G_INVERSE_DECODING);
	dgd_->push_back(std::move(KronDecompFac));

	// one factor of the Kronecker decomposition of G_m^p for arbitrary m.
	KronDecompFac =
		std::make_unique<TransformationZZ_PrimePower>(m, p, TransformationType::INT_G_POWERFUL);
	mgp_->push_back(std::move(KronDecompFac));

	// one factor of the Kronecker decomposition of (G_m^p)^-1 for arbitrary m.
	KronDecompFac =
		std::make_unique<TransformationZZ_PrimePower>(m, p, TransformationType::INT_G_INVERSE_POWERFUL);
	dgp_->push_back(std::move(KronDecompFac));	
}

// Pre-compute all complex vector transformation matrices.
void LweField::computeComplexTransMatrices(pos_int p, pos_int m)
{
	using namespace RLWE_Toolkit::Transformations;

	// one factor of the Kronecker decomposition of CRT_m for arbitrary m.
	std::unique_ptr<TransformationCompCRT::base_type> KronDecompFac =
		std::make_unique<TransformationCompCRT>(m, p, TransformationCompCRT::TransformationType::COMPLEX_CRT_M);
	crtM_->push_back(std::move(KronDecompFac));

	// one factor of the Kronecker decomposition of CRT_m^* for arbitrary m.
	KronDecompFac =
		std::make_unique<TransformationCompCRT>(m, p, TransformationCompCRT::TransformationType::COMPLEX_CRT_M_STAR);
	crtMStar_->push_back(std::move(KronDecompFac));

	// one factor of the Kronecker decomposition of CRT_m^*^-1 for arbitrary m.
	KronDecompFac =
		std::make_unique<TransformationCompCRT>(m, p, TransformationCompCRT::TransformationType::COMPLEX_CRT_M_STAR_INVERSE);
	crtMStarInverse_->push_back(std::move(KronDecompFac));

	// one factor of the Kronecker decomposition of CRT_m^-1 for arbitrary m.
	KronDecompFac =
		std::make_unique<TransformationCompCRT>(m, p, TransformationCompCRT::TransformationType::COMPLEX_CRT_M_INVERSE);
	crtMInverse_->push_back(std::move(KronDecompFac));
}

void LweField::computeRealTransMatrices(pos_int p, pos_int m)
{
	using namespace RLWE_Toolkit::Transformations;

	// one factor of the Kronecker decomposition of D for arbitrary m.
	std::unique_ptr<TransformationRealPrimePower::base_type> KronDecompFac =
		std::make_unique<TransformationRealPrimePower>(m, p, TransformationRealPrimePower::TransformationType::REAL_GS_D);
	gs_decomp_D_->push_back(std::move(KronDecompFac));

	// one factor of the Kronecker decomposition of U for arbitrary m.
	KronDecompFac =
		std::make_unique<TransformationRealPrimePower>(m, p, TransformationRealPrimePower::TransformationType::REAL_GS_U);
	gs_decomp_U_->push_back(std::move(KronDecompFac));

	// one factor of the Kronecker decomposition of C*B' (D) for arbitrary m.
	KronDecompFac =
		std::make_unique<TransformationRealPrimePower>(m, p, TransformationRealPrimePower::TransformationType::REAL_CONVERT_GAUSSIANS);
	sample_gauss_D_->push_back(std::move(KronDecompFac));
}

// Compute the coordinate vector of the neutral element in K for the different considered bases.
void LweField::computeOnes()
{
	pos_int n = getDimension();
	ones = std::make_unique<std::map<RingElement::Basis, std::unique_ptr<coordinate_vec> > >();
	coordinate_vec coords = coordinate_vec(n, 0);
	coords[0] = 1;

	// The first basis element in p is alway 1. So the coordinate vector of 1 in R is (1, 0, ..., 0).
	(*ones)[Basis::BASIS_POWERFUL] = std::make_unique<coordinate_vec>(coords);

	// Computes the 1 in R_q.
	this->applyIntTransformation(TransMat::INTEGER_CRT_MQ, coords, coords);
	(*ones)[Basis::BASIS_CRT] = std::make_unique<coordinate_vec>(coords);

	// The coordinate vector of 1 in H is (1, 1, ..., 1)
	comp_vec vec = comp_vec(n, complex_type(REAL_TYPE_ONE));

	// Compute 1 in R_dual via sigma^-1 (1) in basis d.
	this->applyCompTransformation(TransMat::COMPLEX_CRT_M_STAR, vec, vec);
	convertCompVecToIntVec(coords, vec);
	(*ones)[Basis::BASIS_DECODING] = std::make_unique<coordinate_vec>(coords);

	// convert to t^-1 p
	this->applyIntTransformation(TransMat::INTEGER_L_M_INVERSE, coords, coords);
	(*ones)[Basis::BASIS_T_INVERSE_POWERFUL] = std::make_unique<coordinate_vec>(coords);

	// convert to t^-1 c
	this->applyIntTransformation(TransMat::INTEGER_CRT_MQ, coords, coords);
	(*ones)[Basis::BASIS_T_INVERSE_CRT] = std::make_unique<coordinate_vec>(coords);
}

// Computes the coordinates of the element g in the CRT basis.
void LweField::computeGinCRT(std::list<std::pair<pos_int, pos_int>> primeFactors)
{
	typedef std::list<std::pair<pos_int, pos_int>> PrimeList;
	std::unique_ptr<coordinate_vec> g_ptr = std::make_unique<coordinate_vec>(*(*ones)[Basis::BASIS_CRT]);
	coordinate_vec one = *(*ones)[Basis::BASIS_CRT];

	pos_int pos = 1;
	int q = rlwe_modulus_;

	// traverse primePowers in reversed order.
	for (PrimeList::const_reverse_iterator primePower = primeFactors.crbegin(); primePower != primeFactors.crend(); primePower++)
	{
		pos_int p = primePower->first;
		pos_int m_prime = pow(primePower->first, primePower->second - 1); // p^k-1
		pos_int phi_m = (p - 1) * m_prime; // phi(p)*m' = phi(m)
		if (p == 2){
			continue;
		}
		// coordinates of XI_p in the powerful basis
		coordinate_vec XI_p = coordinate_vec(n_);
		XI_p[m_prime*pos] = 1;

		// convert to CRT basis
		this->applyIntTransformation(TransMat::INTEGER_CRT_MQ, XI_p, XI_p);
		
		// update g
		for (pos_int j = 0; j < n_; j++)
		{
			(*g_ptr)[j] = (*g_ptr)[j] * (one[j] - XI_p[j]);
		}
		*g_ptr = *g_ptr % q;

		pos *= phi_m;
	}

	g_crt_coords_.release();
	g_crt_coords_ = std::move(g_ptr);
}

// Compute the coordinate of t^-1 in the decoding basis.
void LweField::computeT_inverseInD()
{
	// (1, 0, ..., 0) are coordinates of 1 in R in the powerful basis. Thus, (1, 0, ..., 0) are the coordinates
	// of t^-1 in the powerful basis t^-1 p of R^dual.
	std::unique_ptr<coordinate_vec> t_ptr = std::make_unique<coordinate_vec>(*(*ones)[Basis::BASIS_POWERFUL]);

	// convert to decoding basis.
	this->applyIntTransformation(TransMat::INTEGER_L_M_INVERSE, *t_ptr, *t_ptr);
	t_inverse_decoding_coords_.release();
	t_inverse_decoding_coords_ = std::move(t_ptr);
}
