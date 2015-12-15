#pragma once

#ifndef RLC_FIELD_H

#define RLC_FIELD_H

#include <vector>
#include <list>
#include <complex>
#include <tuple>
#include <map>

#include "RLWE/RingLweCryptographyElement.h"
#include "Transformations/AbstractVectorTransformation.h"
#include "Transformations/Transformation algorithms.h"

// typedefs 
typedef int coordinate_type;	// type for coordinate vectors (int or long)
typedef unsigned int pos_int;

typedef double real_type; // type for representation of the reals. Also for usage in std::complex. 
typedef std::complex<real_type> complex_type;

// vecs and matrices
typedef std::vector<coordinate_type> coordinate_vec;
typedef std::vector<real_type> real_vec;
typedef std::vector<complex_type> comp_vec;

class TransformationCompCRT;
class MatrixCompMult;

/*
RingLweCryptographyField.h

This class represents the cyclotomic number field K = QQ(zeta_m), where zeta_m is a primitive m-th root of unity. 
The elements represented by RingLweCryptographyElement have a "live-in" relation to this class. The main tasks 
are the construction of special elements and the management and handling of the needed transformation matrices.
Also stores the LWE modulus q.

This class is based on the descriptions of the companion work:
[May15] Christoph Mayer, Implementing a Toolkit for Ring-LWE Based Cryptography in Arbitrary
		Cyclotomic Number Fields, 2015

@author Christoph Mayer
@version 1.0
*/
class RingLweCryptographyField : public std::enable_shared_from_this<RingLweCryptographyField>
{
public:
	// The vector transformations this class handles.
	enum class TransformationMatrices
	{
		COMPLEX_CRT_M,
		COMPLEX_CRT_M_INVERSE,
		COMPLEX_CRT_M_STAR,
		COMPLEX_CRT_M_STAR_INVERSE,
		INTEGER_L_M,
		INTEGER_L_M_INVERSE,
		INTEGER_CRT_MQ,
		INTEGER_CRT_MQ_INVERSE,
		INTEGER_MULT_G_DEC,
		INTEGER_DIV_G_DEC,
		INTEGER_MULT_G_POW,
		INTEGER_DIV_G_POW,
		REAL_SAMPLE_GAUSS_D,
		REAL_GS_DECOMP_U,
		REAL_GS_DECOMP_D
	};

	/*
	Constructor of K = QQ(zeta_m), where zeta_m is a primitive m-th root of unity. 
	The ring-LWE modulus is given by @param modulus.
	 
	@param m: Represents m-th root of unity.
	@param modulus: r-LWE modulus.
	*/
	RingLweCryptographyField(pos_int m, pos_int modulus);

	// Default destructor
	~RingLweCryptographyField();

	
	
	/**********************\
			Getter
	\**********************/

	inline pos_int getDimension() const { return n_; };
	inline pos_int getModulus() const { return rlwe_modulus_; };
	inline pos_int getM() const { return m_; };
	inline pos_int getRadM() const { return rad_m_; };

	/**********************\
		Special elements
	\**********************/

	// Returns the special element t^-1 in the decoding basis.
	RingLweCryptographyElement getElementT_inverse() const {
		RingElement res = RingElement(shared_from_this(), RingElement::Basis::BASIS_DECODING, *t_inverse_decoding_coords_);
		return res;
	};

	// Returns the special element g in the CRT basis, so g in R_q.
	RingLweCryptographyElement getElementG() const {
		RingElement res = RingElement(shared_from_this(), RingElement::Basis::BASIS_CRT, *g_crt_coords_, rlwe_modulus_);
		return res;
	};

	// Get the zero element in the specific basis.
	RingLweCryptographyElement getZero(RingLweCryptographyElement::Basis basis) const;

	// Get the neutral element in the specific basis.
	RingLweCryptographyElement getOne(RingLweCryptographyElement::Basis basis) const;

	/**********************\
	   Transformations
	\**********************/

	// x = transformation(a) for x,a in ZZ^n
	void applyIntTransformation(TransformationMatrices transformation, coordinate_vec& x, coordinate_vec const& a) const;

	// x = transformation(a) for x,a in RR^n
	void applyRealTransformation(TransformationMatrices transformation, real_vec& x, real_vec const& a) const;

	// x = transformation(a) for x,a in CC^n
	void applyCompTransformation(TransformationMatrices transformation, comp_vec& x, comp_vec const& a) const;
	
private:
	friend class RingLweCryptographyElement;

	pos_int const m_;	// represents m-th cyclotomic root
	pos_int const n_;	// = phi(m_) (Euler totient)
	pos_int rad_m_; // = rad(m) = product of all primes dividing m.

	pos_int const rlwe_modulus_; // ring-LWE modulus for cryptosystems
	
	std::unique_ptr<coordinate_vec> t_inverse_decoding_coords_; // Precomputed coordinates for t^-1 in the decoding basis.
	std::unique_ptr<coordinate_vec> g_crt_coords_; // Precomputed coordinates for g in the CRT basis.

	// Map for precomputed coordinate vectors for the neutral element in different bases.
	std::unique_ptr< std::map<RingElement::Basis, std::unique_ptr<coordinate_vec>> > ones; 

	// unique pointer to vector transformations
	template<typename T> using avt_uptr = std::unique_ptr < AbstractVectorTransformation<T> >;

	// unique pointer to list of unique pointers to vector transformations
	template<typename T> using avt_list_uptr = std::unique_ptr < std::list<avt_uptr<T>> >; 
	
	// Integer transformations
	avt_list_uptr<coordinate_type> Lm_; // Kronecker decomposition of L_m.
	avt_list_uptr<coordinate_type> LmInverse_; // Kronecker decomposition of L_m^-1.

	avt_list_uptr<coordinate_type> crtMq_; // Kronecker decomposition of CRT_m,q.
	avt_list_uptr<coordinate_type> crtMqInverse_; // Kronecker decomposition of CRT_m,q^-1.

	avt_list_uptr<coordinate_type> mgd_; // Kronecker decomposition of matrix for multiplication with g in decoding basis.
	avt_list_uptr<coordinate_type> dgd_; // Kronecker decomposition of matrix for multiplication with g^-1 in decoding basis.
	avt_list_uptr<coordinate_type> mgp_; // Kronecker decomposition of matrix for multiplication with g in powerful basis.
	avt_list_uptr<coordinate_type> dgp_; // Kronecker decomposition of matrix for multiplication with g^-1 in powerful basis.

	// Complex transformations
	avt_list_uptr<complex_type> crtM_; // Kronecker decomposition of CRT_m.
	avt_list_uptr<complex_type> crtMInverse_; // Kronecker decomposition of CRT_m^-1.
	avt_list_uptr<complex_type> crtMStar_; // Kronecker decomposition of CRT_m^*.
	avt_list_uptr<complex_type> crtMStarInverse_; // Kronecker decomposition of CRT_m^*^-1.

	// Real transformations
	avt_list_uptr<real_type> sample_gauss_D_; // Kronecker decomposition of matrix for transforming Gaussians in RR^n into Gaussians in K_RR.
	avt_list_uptr<real_type> gs_decomp_U_; // Kronecker decomposition of matrix U in Gram-Schmidt decomposition of CRT_m.
	avt_list_uptr<real_type> gs_decomp_D_; // Kronecker decomposition of matrix D in Gram-Schmidt decomposition of CRT_m.

	// Init function for construction.
	void init();

	// Pre-compute integral vector transformation matrices.
	// Computes one part of the Kronecker decomposition, the one for the prime power m.
	// m = p^k, m_prime = p^(k-1), n = phi(m).
	void computeIntegerTransMatrices(pos_int p, pos_int m_prime, pos_int m, pos_int n);

	// Pre-compute complex vector transformation matrices.
	// Computes one part of the Kronecker decomposition, the one for the prime power m.
	// m = p^k.
	void computeComplexTransMatrices(pos_int p, pos_int m);

	// Pre-compute real vector transformation matrices.
	// Computes one part of the Kronecker decomposition, the one for the prime power m.
	// m = p^k.
	void computeRealTransMatrices(pos_int p, pos_int m);
	
	// Pre-computes the coordinate vectors for the neutral element in different bases.
	void computeOnes();

	// Pre-computes the coordinate vectors for the special element g in the CRT basis. 
	void computeGinCRT(std::list<std::pair<pos_int, pos_int>> primeFactors);

	// Pre-computes the coordinate vectors for the special element t^-1 in the decoding basis. 
	void computeT_inverseInD();
};

/**********************\
	  Comparison
\**********************/

// Returns true, if and only if, the dimension and the modulus are equal.
inline bool operator==(RingLweCryptographyField const& lhs, RingLweCryptographyField const& rhs) 
{
	return ((lhs.getDimension() == rhs.getDimension()) && (lhs.getModulus() == rhs.getModulus()));
}
inline bool operator!=(RingLweCryptographyField const& lhs, RingLweCryptographyField const& rhs)
{
	return !(lhs == rhs);
}


#endif // !RLC_FIELD_H