#pragma once

#ifndef RLC_ELMT_H

#define RLC_ELMT_H

#include <vector>
#include <complex>
#include <memory>
#include <iostream>

// typedefs 
typedef int coordinate_type;	// type for coordinate vectors (int or long)
typedef unsigned int pos_int;

typedef double real_type; // type for representation of the reals. Also for usage in std::complex. 
typedef std::complex<real_type> complex_type;

// vecs and matrices
typedef std::vector<coordinate_type> coordinate_vec;
typedef std::vector<real_type> real_vec;
typedef std::vector<complex_type> comp_vec;

class RingLweCryptographyField; // forward declaration

/*
RingLweCryptographyElement.h

This class represents elements in the ring of integers R and its dual ideal R^dual in some
cyclotomic number field K. The elements are represented by a coordinate vector with respect
to a specific basis. Furthermore, this class provides a variety of algorithms specialized 
for the usage in ring-LWE based cryptography. 

This class is based on the descriptions of the companion work:
[May15] Christoph Mayer, Implementing a Toolkit for Ring-LWE Based Cryptography in Arbitrary 
		Cyclotomic Number Fields, 2015

@author Christoph Mayer
@version 1.0
*/
class RingLweCryptographyElement
{
public:
	typedef std::shared_ptr<RingLweCryptographyField const> rlwe_field_sptr;

	// Flags for different bases we need.
	enum class Basis
	{
		BASIS_POWERFUL,
		BASIS_T_INVERSE_POWERFUL,
		BASIS_CRT,
		BASIS_T_INVERSE_CRT,
		BASIS_DECODING
	};


	/*
	Default Constructor. Initializes an "empty" element in the field @param field.
	This constructor should only be used when combined with another constructive method, i.e., to produce a new "blank" 
	instance which gets overwritten by a method like RingLweCryptographyField::add or RingLweCryptographyField::mult.
	*/
	RingLweCryptographyElement(rlwe_field_sptr field);

	/*
	Constructor of a known element a in I = R, R^dual. The ideal I in K is represented by the given basis.
	An element in I needs a field K to live in and some coordinates in the given basis. If the
	element a lives actually in the quotient I = I/qI, then q is given as primeModulusQ. It is also 
	possible that the element a lives in a power I^k. Then k is given as idealPower.
	 
	 
	@param field Shared pointer to the field in which the constructed element lives.
	@param basis Represents the current basis.
	@param coordinates The integral coordinate vector of this element. The coordinates are w.r.t. the basis
		   represented by the enum Basis.
	@param primeModulusQ Indicates whether the computation is modulo q. If primeModulusQ > 0, then we set q = primeModulusQ.
	@param idealPower The integer k, if this element is in (R^dual)^k. If the element is in R or R_q, idealPower is set to 0.
	*/
	RingLweCryptographyElement(rlwe_field_sptr field, Basis basis,
		coordinate_vec& coordinates, pos_int primeModulusQ = 0, pos_int idealPower = 1);
	
	// Copy constructor
	RingLweCryptographyElement(RingLweCryptographyElement const& elmt);

	// Default destructor 
	~RingLweCryptographyElement();

	

	/*
	Changes the basis of this element to the new basis and computes the coordinates accordingly.
	
	@param newBasis The new basis in which the element will be represented in. 
	*/
	bool changeBasisTo(Basis newBasis);


	/************************************************************************\

							Operator overloading

	\************************************************************************/

	RingLweCryptographyElement& operator=(RingLweCryptographyElement const& elmt);

	RingLweCryptographyElement& operator+=(RingLweCryptographyElement const& elmt);
	RingLweCryptographyElement& operator-=(RingLweCryptographyElement const& elmt);

	RingLweCryptographyElement& operator*=(RingLweCryptographyElement const& elmt);
	RingLweCryptographyElement& operator*=(const int scalar);


	/********************************************************************************\
	
									Getter & Setter

	\********************************************************************************/

	// Returns a reference to the actual field object.
	inline RingLweCryptographyField const& getField() const { return *field_; }; 

	// Returns the shared pointer pointing to the field object.
	inline rlwe_field_sptr getFieldSptr() const { return field_; }; 
	
	inline coordinate_vec& getCoordinates() const { return *coordinates_; };
	inline void setCoordinates(coordinate_vec& coordinates){ *coordinates_ = coordinates; };

	inline Basis getBasis() const { return basis_; };
	inline void setBasis(Basis basis){ basis_ = basis; };

	inline pos_int getPrimeModulusQ() const { return primeModulusQ_; };
	inline void setPrimeModulusQ(pos_int q){ primeModulusQ_ = q; };

	inline pos_int getIdealPower() const { return k_; };
	inline void setIdealPower(pos_int idealPower){
		if (basis_ == Basis::BASIS_POWERFUL || basis_ == Basis::BASIS_CRT){
			k_ = 0;
		}
		else{
			k_ = idealPower;
		}
	};


private:

	/*********************\
		class variables
	\*********************/

	std::shared_ptr<RingLweCryptographyField const> field_; // The field in which this element lives.

	std::unique_ptr<coordinate_vec> coordinates_;		// Integral coordinate vector of the element w.r.t. basis.

	Basis basis_;						// The actual basis in which the element is represented.

	pos_int primeModulusQ_;		// Possible prime modulus q. If primeModulusQ > 0, then q = primeModulusQ
										// and all computations are mod q.

	pos_int k_;		// Represents I = (R^dual)^k, default k_ = 1 or k_ = 0 if this element is in R or R_q.

	/*********************************\
				Friends
	\*********************************/

	friend bool add(RingLweCryptographyElement& x, RingLweCryptographyElement const& a, RingLweCryptographyElement const& b);
	friend bool mult(RingLweCryptographyElement& x, RingLweCryptographyElement const& a, RingLweCryptographyElement const& b);

	/*********************************\
		Help functions addition
	\*********************************/

	// Adds the entries of a and b component-wise and stores the result in x.
	template<typename T> static void addComponentWise(std::vector<T>& x, std::vector<T> const& a, std::vector<T> const& b);
	// Adds the entries of a and b component-wise modulo modulus and stores the result in x.
	static void addComponentWiseModulo(coordinate_vec& x, coordinate_vec const& a, coordinate_vec const& b, int modulus);
	
	/*
	Adds the elements a and b via the canonical embedding sigma. That is, the elements are embedded into H, added component-wise 
	and then pulled back into K. This is only necessary if one element is in R and the other in R^dual. The result is stored in x.
	As an intermediate step we have to compute with vectors of complex numbers. These are then converted into integral vectors, 
	as the stored complex values are very close to an integer. The conversion actually rounds off the real part of the complex number.
	*/
	static void addViaEmbedding(RingLweCryptographyElement& x, RingLweCryptographyElement const& a, 
		RingLweCryptographyElement const& b);
	
	// Subroutine for addition in the case primeModulusQ is greater 0 and addition is performed modulo.
	static bool addModulo(RingLweCryptographyElement& x, RingLweCryptographyElement& a, RingLweCryptographyElement& b);
	// Subroutine for addition in the case primeModulusQ is 0 and addition is performed non modulo.
	static bool addNonModulo(RingLweCryptographyElement& x, RingLweCryptographyElement& a, RingLweCryptographyElement& b);

	// Small subroutine. Adds the elements a and b and sets x.basis_ to basis.
	static void addAndSetBasis(RingLweCryptographyElement& x, RingLweCryptographyElement const& a, 
		RingLweCryptographyElement const& b, Basis basis);

	/***************************************\
		Help functions multiplication
	\***************************************/
	
	// Multiplies the entries of a and b component-wise and stores the result in x.
	template<typename T> static void multComponentWise(std::vector<T>& x, std::vector<T> const& a, std::vector<T> const& b);
	// Multiplies the entries of a and b component-wise modulo modulus and stores the result in x.
	static void multComponentWiseModulo(coordinate_vec& x, coordinate_vec const& a, coordinate_vec const& b, int modulus);

	// Subroutine for multiplication in the case primeModulusQ is greater 0 and multiplication is performed modulo.
	static bool multModulo(RingLweCryptographyElement& x, RingLweCryptographyElement& a, RingLweCryptographyElement& b);
	// Subroutine for multiplication in the case primeModulusQ is 0 and multiplication is performed non modulo.
	static bool multNonModulo(RingLweCryptographyElement& x, RingLweCryptographyElement& a, RingLweCryptographyElement& b);

	/*
	Multiplies the elements a and b via the canonical embedding sigma. That is, the elements are embedded into H, multiplied component-wise 
	and then pulled back into K. The result is stored in x. As an intermediate step we have to compute with vectors of complex numbers.
	These are then converted into integral vectors, as the stored complex values are very close to an integer. The conversion actually
	rounds off the real part of the complex number.
	*/
	static void multiplyViaEmbedding(RingLweCryptographyElement& x, RingLweCryptographyElement const& a, RingLweCryptographyElement const& b);
	// Multiplies the coordinate vectors component-wise.
	static void multiplyInCrt(RingLweCryptographyElement& x, RingLweCryptographyElement const& a, RingLweCryptographyElement const& b);


	/***************************************\
		Help functions for base switching
		& output basis computation
	\***************************************/
	
	// Subroutines for the change basis operation. Changes the basis by multiplying the coordinate vector with the
	// right transformation matrix. Returns true if successful, false otherwise. 

	// @Param direction: true: p -> c
	//					 false:: c -> p
	bool swapBasesPowerfulAndCrt(bool direction);
	
	// @Param direction: true: t^-1*p -> d
	//					 false:: d -> t^-1*p
	bool swapBasesTInversePowerfulAndDecoding(bool direction);	
};

typedef RingLweCryptographyElement RingElement;

/******************************************************************************\

									Comparison

\******************************************************************************/

// Returns true, if and only if, all class variables of lhs and 
// rhs return true via the respectively == operator.
bool operator==(RingElement const& lhs, RingElement const& rhs);
bool operator!=(RingElement const& lhs, RingElement const& rhs);

/******************************************************************************\

									Addition

\******************************************************************************/	

// Operator overloading:

RingElement operator+(RingElement const& a, RingElement const& b); // Addition
RingElement operator-(RingElement const& a, RingElement const& b); // Subtraction
RingElement operator-(RingElement const& a); // unary minus

// Procedural versions:

/*
Adds two elements a and b and stores the result in x. Thus, x = a + b.
This design leaves the decision, whether to produce a new instance or reuse others, to the user, as
a call like add(a,a,b) is also possible. New "empty" instances can be produces by the 
constructor RingLweCryptographyElement(rlwe_field_sptr field).

@param x The element where the result is stored.
@param a First element in the addition.
@param b Second element in the addition.
*/
bool add(RingElement& x, RingElement const& a, RingElement const& b);

/*
Subtracts the element b from a and stores the result in x. Thus x = a - b.
This design leaves the decision, whether to produce a new instance or reuse others, to the user, as
a call like sub(a,a,b) is also possible. New "empty" instances can be produces by the
constructor RingLweCryptographyElement(rlwe_field_sptr field).

@param x The element where the result is stored.
@param a First element in the subtraction.
@param b Second element in the subtraction.
*/
bool sub(RingElement& x, RingElement const& a, RingElement const& b);

/*
Stores -a in x. Thus x = -a.
This design leaves the decision, whether to produce a new instance or reuse others, to the user, as
a call like neg(a,a) is also possible. New "empty" instances can be produces by the
constructor RingLweCryptographyElement(rlwe_field_sptr field).

@param x The element where the result is stored.
@param a Element to be negated.
*/
bool neg(RingElement& x, RingElement const& a);


/******************************************************************************\

							   Multiplication

\******************************************************************************/

// Operator overloading:

RingElement operator*(RingElement const& a, RingElement const& b); // binary multiplication
RingElement operator*(coordinate_type const scalar, RingElement const& b); // scalar multiplication

// Procedural versions:

/*
Multiplies two elements a and b and stores the result in x. Thus x = a * b.
This design leaves the decision, whether to produce a new instance or reuse others, to the user, as
a call like mult(a,a,b) is also possible. New "empty" instances can be produces by the 
constructor RingLweCryptographyElement(rlwe_field_sptr field).

@param x The element where the result is stored.
@param a First element in the Multiplication.
@param b Second element in the Multiplication.
*/
bool mult(RingElement& x, RingElement const& a, RingElement const& b);

/*
Multiplies the elements a with a scalar and stores the result in x. Thus x = scalar * a.
This design leaves the decision, whether to produce a new instance or reuse others, to the user, as
a call like mult(a,scalar,a) is also possible. New "empty" instances can be produces by the 
constructor RingLweCryptographyElement(rlwe_field_sptr field).

@param x The element where the result is stored.
@param scalar The integer scalar a is multiplied by.
@param a The ring element.
*/
bool multWithScalar(RingElement& x, coordinate_type const scalar, RingElement const& a);


/******************************************************************************\

							Special multiplication

\******************************************************************************/

/*
Computes x = t * a. As R^dual = <t^-1>, multiplying an element a in R^dual with t leads to an element in R.
This is not covered by the normal multiplication routine and is thus separately implemented. Returns true,
if a is an element in R^dual, false otherwise. Sores the result in x.

@param x The result is stored in x. Thus x = t * a.
@param a The element that is multiplied with t.

@return true if successful, false otherwise.
*/
bool multWithT(RingElement& x, RingElement const& a);

/*
Computes x = t^-1 * a. As R^dual = <t^-1>, multiplying an element a in R with t^-1 leads to an element in R^dual.
This can be done by actually multiplying a with the element t^-1 which is precomputed at the construction of a 
RingLweCryptographyField or by switching the bases in the right manner. The second attempt is more efficient
as it only has to deal with integer arithmetic, but not always possible. This routine decides which way is the
best to use. 

@param x The result is stored in x. Thus x = t^-1 * a.
@param a The element that is multiplied with t^-1.

@return true if successful, false otherwise.
*/
bool multWithT_Inverse(RingElement& x, RingElement const& a);


/******************************************************************************\

								Modulo operations

\******************************************************************************/

// Procedural versions:

/*
Computes x = a mod q. That is the coordinate vector of a is taken component-wise modulo q.
Thereby 0 <= z mod q < q for all integers z.

@param x The result is stored in x. Thus x = a mod q.
@param a The element that is taken mod q.
@param q The modulus.
*/
void mod(RingElement& x, RingElement const& a, int q);

/*
Computes x = a mod q. That is the vector is taken component-wise modulo q.
Thereby 0 <= z mod q < q for all integers z.

@param x The result is stored in x. Thus x = a mod q.
@param a The vector that is taken component-wise mod q.
@param q The modulus.
*/
void mod(coordinate_vec& x, coordinate_vec const& a, int q);

// Operator versions:

RingElement operator%(RingElement const& a, int q); // modulo for ring elements
coordinate_vec operator%(coordinate_vec const& a, int q); // modulo for integral vectors


/******************************************************************************\

					Sampling, Discretizing, Decoding

\******************************************************************************/

/*
Decodes the element a into x as described in the companion work [May15]. 
A call of decode(x,a) makes only sense if the element a is in R^dual_q.
*/
void decode(RingElement& x, RingElement const& a);

/*
If a has an coordinate vector in ZZq, i.e., is in some ideal R_q, R^dual_q,..., the procedure 
stores the congruent element with coordinates in [-q/2, q/2] in x and returns true, otherwise 
it returns false. 

@param x The unique representative of a is stored in x.
@param a An element in R_q, R^dual_q,... .
*/
bool computeUniqueRepresentative(RingElement& x, RingElement const& a);

/*
Let b be the basis of the coset representative. Then this routine discretizes an arbitrary element r in K_R, given 
by realCoordinateVector in the basis b to a coset c + p * L(B) of the lattice L(B). Thereby, the lattice L(B) is 
the embedding of the ideal in K represented by b, thus L(B) = sigma(<b>) or equivalent B = sigma(b). The scalar p is 
given by scalingFactor.

@param x The discretized element with coordinate vector w.r.t. basis of cosetRepresentative is stored in x.
@param coordinateVector The coordinate Vector of the element which is discretized.
@param cosetRepresentative The representative of the lattice coset on which the element x is discretized to.
@param scalingFactor A scaling factor for the lattice. Default is 1.
*/
void discretize(RingElement& x, real_vec& realCoordinateVector, RingElement const& cosetRepresentative, 
	int const scalingFactor = 1);

/*
Samples a coordinate vector for an element x in K_RR w.r.t. the decoding basis. The sampled elements
are Gaussian distributed in K_RR.

@param x The coordinate vector of a Gaussian in K_RR w.r.t. the decoding basis.
@param field The field K.
@param mean Mean for the Gaussian distribution.
@param stddev Standard deviation for the Gaussian distribution.
*/
void sampleGaussiansInK_RR(real_vec& x, RingLweCryptographyField const& field, real_type mean, real_type stddev);

/*
Sampling from the discretized Gaussian distribution over c + p*R_dual, 
where c is the coset representative p = scalingFactor. Since the routine 
samples Gaussians in K_RR via sampleGaussiansInK_RR, which is done in 
the decoding basis, c needs to be represented also in the decoding basis.
The computed sample from K_RR is then discretized to c + p*R.

@param x The sampled element is stored in x.
@param mean Mean for the Gaussian distribution.
@param stddev Standard deviation for the Gaussian distribution.
@param cosetRepresentative Coset representative for the discretization.
@param scalingFactor Scaling factor for the discretization. Default is 1.
*/
void sampleDiscretizedGaussian(RingElement& x, real_type mean, real_type stddev, RingElement const& cosetRepresentative,
	int const scalingFactor = 1);

/*
Sampling from the discrete Gaussian distribution over ZZ via rejection sampling.
This algorithm follows the work of Gentry and Peikert. For detailed explanation and analysis I refer to

How to Use a Short Basis: Trapdoors for Hard Lattices and New Cryptographic Constructions,
C. Peikert, C. Gentry, 2008.

@param stddev Standard deviation for the Gaussian distribution.
@param center Center / mean for the Gaussian distribution.
@param funcValue Elements are sampled in [center - stddev*funcValue, center + stddev*funcValue]. Default is 3 
*	(99.73% of the values of the Gaussian distribution are within 3 times the standard deviation from the center).
*/
coordinate_type sampleDiscreteGaussianInZZ(real_type stddev, real_type center, real_type funcValue = 3);

/*
Sampling from the discrete Gaussian distribution over R + c, a coset of the ring of algebraic integers. 
The mean is always 0. This algorithm follows the work of Gentry and Peikert. For detailed explanation 
and analysis I refer to

How to Use a Short Basis: Trapdoors for Hard Lattices and New Cryptographic Constructions,
C. Peikert, C. Gentry, 2008.

@param x Sampled element is stored in x.
@param c Representative of the coset c + R.
@param stddev Standard deviation for the Gaussian distribution.
*/
void sampleDiscreteGaussianInR(RingElement& x, RingElement const& c, real_type stddev);

/*
Samples a uniformly distributed coordinate vector for the element x. Only the coordinates of x get changed.
The entries are uniformly distributed over [lb, ub].

@param x The element a new coordinate vector is sampled for.
@param lb Lower bound for Interval where the entries are sampled from.
@param ub Upper bound for Interval where the entries are sampled from.
*/
void sampleUniformlyCoordVec(RingElement& x, coordinate_type lb, coordinate_type ub);

/***************************\
		  Output
\***************************/

// Print operator for ring elements.
std::ostream& operator<<(std::ostream& out, RingLweCryptographyElement const& elmt);

/********************************** Help-functions to switch between ZZ^n and CC^n **********************************/

// Converts a vector of integers into a vector of complex numbers.
void convertIntVecToCompVec(comp_vec& compVec, coordinate_vec const& intVec);

// Converts a vector of complex numbers into a vector of integers. We assume that all the complex numbers are very close
// to an integer. The errors are eliminated by the round procedure.
void convertCompVecToIntVec(coordinate_vec& intVec, comp_vec const& compVec);


/********************************** Implementation of template functions ********************************************/

// Component-wise addition of two vectors a, b of type T. Stores the result in x.
template <typename T> void RingElement::addComponentWise(std::vector<T>& x, std::vector<T> const& a, std::vector<T> const& b)
{
	std::transform(a.cbegin(), a.cend(), b.cbegin(), x.begin(), std::plus<T>());
}

// Component-wise multiplication of two vectors a, b of type T. Stores the result in x.
template <typename T> void RingElement::multComponentWise(std::vector<T>& x, std::vector<T> const& a, std::vector<T> const& b)
{
	std::transform(a.cbegin(), a.cend(), b.cbegin(), x.begin(), std::multiplies<T>());
}


#endif // !RLC_ELMT_H