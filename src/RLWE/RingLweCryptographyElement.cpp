#define _USE_MATH_DEFINES


#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <numeric>

#include "RLWE/RingLweCryptographyElement.h"
#include "RLWE/RingLweCryptographyField.h"

#define REAL_TYPE_ZERO 0.0
#define REAL_TYPE_ONE 1.0
#define COMP_M_I complex_type(REAL_TYPE_ZERO, REAL_TYPE_ONE)
#define COMP_M_SQRT1_2 complex_type(M_SQRT1_2)

typedef RingLweCryptographyElement RingElement;
typedef RingLweCryptographyField LweField;
typedef RingElement::Basis Basis;
typedef RingLweCryptographyField::TransformationMatrices TransMat;

// Default constructor
RingElement::RingLweCryptographyElement(rlwe_field_sptr field) 
	:
	field_(field),
	coordinates_(std::make_unique<coordinate_vec>(field_->getDimension(), 0)),
	basis_(Basis::BASIS_POWERFUL),
	primeModulusQ_(0),
	k_(0)
{}

// Constructs an element in the given field, with coordinates in the given basis 
// and optional values for the prime modulus q and the ideal power k.
RingElement::RingLweCryptographyElement(rlwe_field_sptr field, RingElement::Basis basis,
	coordinate_vec& coordinates, pos_int primeModulusQ, pos_int idealPower)
	:
	field_(field),
	coordinates_(std::make_unique<coordinate_vec>(coordinates)),
	basis_(basis),
	primeModulusQ_(primeModulusQ),
	k_(idealPower)
{
	if (basis_ == Basis::BASIS_POWERFUL || basis_ == Basis::BASIS_CRT){
		// represented ideal is R or R_q and k can only be zero
		k_ = 0;
	}
	if (primeModulusQ_ > 0){ // make sure that coordinate entries are between 0 and q
		*coordinates_ = *coordinates_ % primeModulusQ_;
	}
}

// Copy constructor
RingElement::RingLweCryptographyElement(RingElement const& elmt)
	:
	field_(elmt.getFieldSptr()),
	coordinates_(std::make_unique<coordinate_vec>(elmt.getCoordinates())),
	basis_(elmt.getBasis()),
	primeModulusQ_(elmt.getPrimeModulusQ()),
	k_(elmt.getIdealPower())
{}

// Default destructor
RingElement::~RingElement()
{}

bool operator==(RingElement const& lhs, RingElement const& rhs)
{
	return ((lhs.getField() == rhs.getField())
		&& (lhs.getBasis() == rhs.getBasis())
		&& (lhs.getCoordinates() == rhs.getCoordinates())
		&& (lhs.getPrimeModulusQ() == rhs.getPrimeModulusQ())
		&& (lhs.getIdealPower() == rhs.getIdealPower()));
}

bool operator!=(RingElement const& lhs, RingElement const& rhs)
{
	return !(lhs == rhs);
}

RingElement& RingElement::operator=(RingElement const& elmt)
{
	this->field_ = elmt.getFieldSptr();
	this->coordinates_ = std::make_unique<coordinate_vec>(elmt.getCoordinates());
	this->basis_ = elmt.getBasis();
	this->primeModulusQ_ = elmt.getPrimeModulusQ();
	this->k_ = elmt.getIdealPower();
	return *this;
}

/************************************************************************************************************\
					
											    Addition

\************************************************************************************************************/

// x = a + b
bool add(RingElement& x, RingElement const& a, RingElement const& b)
{
	// possible cases: modulus = max(0,0), modulus = max(q,0), 
	// modulus = max(0,q), modulus = max(q,q)
	pos_int modulus = std::max(a.getPrimeModulusQ(), b.getPrimeModulusQ());

	// addition is always performed in the bigger ideal. For k>=l it is 
	// I^k subset I^l
	pos_int k = std::min(a.getIdealPower(), b.getIdealPower());
	bool isModulo = modulus > 0;

	// create working copy
	RingElement a_copy = a;
	RingElement b_copy = b;

	// check if the elements live in the same field
	if (x.getField() != a.getField()){
		throw "Fields mismatch!";
	}
	x.setPrimeModulusQ(modulus);
	
	if (isModulo){ // addition modulo q
		if (RingElement::addModulo(x, a_copy, b_copy)){
			x.setIdealPower(k);
			return true;
		}
		return false;
	}
	else
	{ // normal addition
		if (RingElement::addNonModulo(x, a_copy, b_copy)){
			x.setIdealPower(k);
			return true;
		}
		return false;
	}
}

// x = a - b
bool sub(RingElement& x, RingElement const& a, RingElement const& b)
{
	return neg(x, b) &&	add(x, a, x);
}

// x = -a
bool neg(RingElement& x, RingElement const& a)
{
	x = (-1) * a;

	int q = x.getPrimeModulusQ();
	if (q > 0){ // negation modulo q
		coordinate_vec& coords = x.getCoordinates();
		std::for_each(coords.begin(), coords.end(), [q](coordinate_type& y){ y += q; });
	}
	
	return true;
}


/***********************************************\
                   Operators
\***********************************************/

RingElement operator+(RingElement const& a, RingElement const& b)
{
	RingElement x(a.getFieldSptr());
	add(x, a, b);
	return x;
}

RingElement operator-(RingElement const& a, RingElement const& b)
{
	return a + (-b);
}

RingElement operator-(RingElement const& a)
{
	RingElement x(a.getFieldSptr());
	neg(x, a);
	return x;
}

RingElement& RingElement::operator+=(RingElement const& elmt)
{
	add(*this, *this, elmt);
	return *this;
}

RingElement& RingElement::operator-=(RingElement const& elmt)
{
	sub(*this, *this, elmt);
	return *this;
}


/***********************************************\
            Help functions addition
\***********************************************/

// Algorithm for addition x = a + b in the modulo case. Looks up the two bases, transforms them accordingly, adds the 
// coordinate vectors mod primeModulusQ and updates the output basis.
bool RingElement::addModulo(RingElement& x, RingElement& a, RingElement& b)
{
	Basis b1 = a.getBasis();
	Basis b2 = b.getBasis();

	if (b1 == b2){
		RingElement::addAndSetBasis(x, a, b, b1);
		return true;
	}

	switch (b1)
	{
	case Basis::BASIS_POWERFUL:
		a.swapBasesPowerfulAndCrt(true);
	case Basis::BASIS_CRT:
		switch (b2)
		{
		case Basis::BASIS_POWERFUL:
			b.swapBasesPowerfulAndCrt(true);
		case Basis::BASIS_CRT:					// (p, c) 
			RingElement::addAndSetBasis(x, a, b, Basis::BASIS_CRT);
			return true;

		case Basis::BASIS_DECODING:				// (p, d) , (c, d)
			b.swapBasesTInversePowerfulAndDecoding(false);
		case Basis::BASIS_T_INVERSE_POWERFUL:	// (p, t^-1 p) , (c, t^-1 p)
			b.swapBasesPowerfulAndCrt(true);
		case Basis::BASIS_T_INVERSE_CRT:		// (p, t^-1 c) , (c, t^-1 c)
			RingElement::addAndSetBasis(x, a, b, Basis::BASIS_T_INVERSE_CRT);
			return true;

		default:
			return false;
		}

	case Basis::BASIS_DECODING:
		a.swapBasesTInversePowerfulAndDecoding(false);
	case Basis::BASIS_T_INVERSE_POWERFUL:
		a.swapBasesPowerfulAndCrt(true);
	case Basis::BASIS_T_INVERSE_CRT:
		switch (b2)
		{
		case Basis::BASIS_POWERFUL:				// (d, p) , (t^-1 p, p) , (t^-1 c, p)
			b.swapBasesPowerfulAndCrt(true);
		case Basis::BASIS_CRT:					// (d, c) , (t^-1 p, c) , (t^-1 c, c)
			RingElement::addAndSetBasis(x, a, b, Basis::BASIS_T_INVERSE_CRT);
			return true;

		case Basis::BASIS_DECODING:				// (t^-1 p, d) , (t^-1 c, d)
			b.swapBasesTInversePowerfulAndDecoding(false);
		case Basis::BASIS_T_INVERSE_POWERFUL:	// (d, t^-1 p) , (t^-1 c, t^-1 p) 
			b.swapBasesPowerfulAndCrt(true);
		case Basis::BASIS_T_INVERSE_CRT:		// (d, t^-1 c) , (t^-1 p, t^-1 c) 	
			RingElement::addAndSetBasis(x, a, b, Basis::BASIS_T_INVERSE_CRT);
			return true;

		default:
			return false;
		}

	default:
		return false;
	}
}

// Algorithm for addition x = a + b in the non-modulo case. Looks up the two bases, transforms them accordingly, adds the 
// coordinate vectors and updates the output basis.
bool RingElement::addNonModulo(RingElement& x, RingElement& a, RingElement& b)
{
	Basis b1 = a.getBasis();
	Basis b2 = b.getBasis();

	if (b1 == b2){
		RingElement::addAndSetBasis(x, a, b, b1);
		return true;
	}

	switch (b1)
	{
	case Basis::BASIS_POWERFUL:
		switch (b2)
		{
		case Basis::BASIS_T_INVERSE_POWERFUL:	// (p, t^-1 p)
			b.swapBasesTInversePowerfulAndDecoding(true);
		case Basis::BASIS_DECODING:				// (p, d)
			RingElement::addViaEmbedding(x, a, b);
			x.setBasis(Basis::BASIS_DECODING);
			return true;

		default:
			return false;
		}

	case Basis::BASIS_T_INVERSE_POWERFUL:
		a.swapBasesTInversePowerfulAndDecoding(true);
	case Basis::BASIS_DECODING:
		switch (b2)
		{
		case Basis::BASIS_POWERFUL:				// (t^-1 p, p) , (d, p)
			RingElement::addViaEmbedding(x, b, a);
			x.setBasis(Basis::BASIS_DECODING);
			return true;

		case Basis::BASIS_T_INVERSE_POWERFUL:	// (d, t^-1 p)	
			b.swapBasesTInversePowerfulAndDecoding(true);
		case Basis::BASIS_DECODING:				// (t^-1 p, d)
			RingElement::addAndSetBasis(x, a, b, Basis::BASIS_DECODING);
			return true;

		default:
			return false;
		}
	}
}

// Adds the two vectors a, b component-wise modulo @param modulus and stores the result in x.
void RingElement::addComponentWiseModulo(coordinate_vec& x, coordinate_vec const& a, coordinate_vec const& b, int modulus)
{
	std::transform(a.cbegin(), a.cend(), b.cbegin(), x.begin(), 
		[modulus](const int& x, const int& y){ return (x + y) % modulus; });
}

// Small help function: adds component-wise and updates the basis of x.
void RingElement::addAndSetBasis(RingElement& x, RingElement const& a, RingElement const& b, Basis basis)
{
	int q = x.getPrimeModulusQ();
	if (q > 0){
		RingElement::addComponentWiseModulo(x.getCoordinates(), a.getCoordinates(), b.getCoordinates(), q);
	}
	else
	{
		RingElement::addComponentWise<coordinate_type>(x.getCoordinates(), a.getCoordinates(), b.getCoordinates());
	}
	x.setBasis(basis);
}

// Adds two elements a, b via the embedding sigma. That is a + b = sigma^-1 ( sigma(a) + sigma(b) ). This is just for the
// special case when a.basis_ == BASIS_POWERFUL and b.basis_ == BASIS_DECODING. These bases are expected without further
// testing.
void RingElement::addViaEmbedding(RingElement& x, RingElement const& a, RingElement const& b)
{
	RingLweCryptographyField const& field = a.getField();
	int n = field.getDimension();
	comp_vec sigmaA(n);
	comp_vec sigmaB(n);

	// convert a and into complex vecs
	comp_vec aComp = comp_vec(n);
	convertIntVecToCompVec(aComp, a.getCoordinates());
	comp_vec bComp = comp_vec(n);
	convertIntVecToCompVec(bComp, b.getCoordinates());

	// embed a and b into H
	field.applyCompTransformation(TransMat::COMPLEX_CRT_M, sigmaA, aComp);
	field.applyCompTransformation(TransMat::COMPLEX_CRT_M_STAR_INVERSE, sigmaB, bComp);

	// component-wise addition in H
 	RingElement::addComponentWise<complex_type>(sigmaA, sigmaA, sigmaB);

	// recover x = a + b via the inverse embedding
	field.applyCompTransformation(TransMat::COMPLEX_CRT_M_STAR, aComp, sigmaA);
	convertCompVecToIntVec(x.getCoordinates(), aComp);
}


/********************************************************************************************\
				
				                        Multiplication
				
\********************************************************************************************/

// x = a * b
bool mult(RingElement& x, RingElement const& a, RingElement const& b)
{
	// possible cases: modulus = max(0,0), modulus = max(q,0), 
	// modulus = max(0,q), modulus = max(q,q)
	pos_int modulus = std::max(a.getPrimeModulusQ(), b.getPrimeModulusQ());
	
	bool isModulo = modulus > 0;

	// create working copy
	RingElement a_copy = a;
	RingElement b_copy = b;

	// check if the elements live in the same field
	if (x.getField() != a.getField()){
		x = RingLweCryptographyElement(a.getFieldSptr());
	}
	x.setPrimeModulusQ(modulus);

	if (isModulo){// multiplication modulo q
		return RingElement::multModulo(x, a_copy, b_copy);
	}
	else
	{// normal multiplication
		return RingElement::multNonModulo(x, a_copy, b_copy);
	}
}

// x = scalar * b
bool multWithScalar(RingElement& x, const coordinate_type scalar, RingElement const& a)
{
	int n = a.getField().getDimension();
	if (x != a){
		if (x.getField() != a.getField()){
			return false;
		}
		if (x.getCoordinates().size() != n){
			x.getCoordinates().resize(n);
		}
		x.setBasis(a.getBasis());
		x.setIdealPower(a.getIdealPower());
	}
	coordinate_vec& coord = a.getCoordinates();
	if (a.getPrimeModulusQ() > 0){// scalar multiplication modulo q
		int q = a.getPrimeModulusQ();
		std::transform(coord.cbegin(), coord.cend(), x.getCoordinates().begin(),
			[q, scalar](const int& entry){ return (entry * scalar) % q; });
		x.setPrimeModulusQ(q);
	}
	else{// multiply component-wise with scalar
		std::transform(coord.cbegin(), coord.cend(), x.getCoordinates().begin(),
			std::bind2nd(std::multiplies<coordinate_type>(), scalar));
	}
	return true;
}


/***********************************************\
                   Operators
\***********************************************/

RingElement operator*(RingElement const& a, RingElement const& b)
{
	RingElement x(a.getFieldSptr());
	mult(x, a, b);
	return x;
}

RingElement operator*(const int scalar, RingElement const& b)
{
	RingElement x(b.getFieldSptr());
	multWithScalar(x, scalar, b);
	return x;
}

RingElement& RingElement::operator*=(RingElement const& elmt)
{
	mult(*this, *this, elmt);
	return *this;
}

RingElement& RingElement::operator*=(const int scalar)
{
	multWithScalar(*this, scalar, *this);
	return *this;
}

/***********************************************\
         Help functions multiplication
\***********************************************/

// Algorithm for multiplication x = a * b in the modulo case. Looks up the two bases, transforms them accordingly, 
// multiplies the coordinate vectors mod primeModulusQ and updates the output basis.
bool RingElement::multModulo(RingElement& x, RingElement& a, RingElement& b)
{
	Basis b1 = a.getBasis();
	Basis b2 = b.getBasis();

	// if a and b are in R, then k = 0 + 0 = 0
	// if a in R and b in R^dual^k', then k = 0 + k' = k'
	// if a in R^dual^s and b in R^dual^s, then k = s + t
	// as desired
	unsigned int k = a.getIdealPower() + b.getIdealPower();

	switch (b1)
	{
	case Basis::BASIS_POWERFUL:
		a.swapBasesPowerfulAndCrt(true);
	case Basis::BASIS_CRT:
		switch (b2)
		{
		case Basis::BASIS_POWERFUL:				// (p, p) , (c, p)
			b.swapBasesPowerfulAndCrt(true);
		case Basis::BASIS_CRT:					// (p, c) , (c, c)
			RingElement::multiplyInCrt(x, a, b);
			x.setBasis(Basis::BASIS_CRT);
			break;

		case Basis::BASIS_DECODING:				// (p, d) , (c, d)
			b.swapBasesTInversePowerfulAndDecoding(false);
		case Basis::BASIS_T_INVERSE_POWERFUL:	// (p, t^-1 p) , (c, t^-1 p)
			b.swapBasesPowerfulAndCrt(true);
		case Basis::BASIS_T_INVERSE_CRT:		// (p, t^-1 c) , (c, t^-1 c)
			RingElement::multiplyInCrt(x, a, b);
			x.setBasis(Basis::BASIS_T_INVERSE_CRT);
			break;

		default:
			return false;
		}
		break;

	case Basis::BASIS_DECODING:
		a.swapBasesTInversePowerfulAndDecoding(false);
	case Basis::BASIS_T_INVERSE_POWERFUL:
		a.swapBasesPowerfulAndCrt(true);
	case Basis::BASIS_T_INVERSE_CRT:
		switch (b2)
		{
		case Basis::BASIS_POWERFUL:				// (d, p) , (t^-1 p, p) , (t^-1 c, p)
			b.swapBasesPowerfulAndCrt(true);
		case Basis::BASIS_CRT:					// (d, c) , (t^-1 p, c) , (t^-1 c, c)
			RingElement::multiplyInCrt(x, a, b);
			x.setBasis(Basis::BASIS_T_INVERSE_CRT);
			break;

		case Basis::BASIS_DECODING:				// (d, d) , (t^-1 p, d) , (t^-1 c, d)
			b.swapBasesTInversePowerfulAndDecoding(false);
		case Basis::BASIS_T_INVERSE_POWERFUL:	// (d, t^-1 p) , (t^-1 p, t^-1 p) , (t^-1 c, t^-1 p)
			b.swapBasesPowerfulAndCrt(true);
		case Basis::BASIS_T_INVERSE_CRT:		// (d, t^-1 c) , (t^-1 p, t^-1 c) , (t^-1 c, t^-1 c)	
			RingElement::multiplyInCrt(x, a, b);
			x.setBasis(Basis::BASIS_T_INVERSE_CRT);
			break;

		default:
			return false;
		}
		break;

	default:
		return false;
	}

	x.setIdealPower(k);
	return true;
}

// Algorithm for multiplication x = a * b in the non-modulo case. Looks up the two bases, transforms them accordingly, 
// multiplies the coordinate vectors and updates the output basis.
bool RingElement::multNonModulo(RingElement& x, RingElement& a, RingElement& b)
{
	Basis b1 = a.getBasis();
	Basis b2 = b.getBasis();

	// if a and b are in R, then k = 0 + 0 = 0
	// if a in R and b in R^dual^k', then k = 0 + k' = k'
	// if a in R^dual^s and b in R^dual^s, then k = s + t
	// as desired
	unsigned int k = a.getIdealPower() + b.getIdealPower();

	switch (b1)
	{
	case Basis::BASIS_POWERFUL:
		switch (b2)
		{
		case Basis::BASIS_POWERFUL:				// (p, p)
			RingElement::multiplyViaEmbedding(x, a, b);
			x.setBasis(Basis::BASIS_POWERFUL);
			break;

		case Basis::BASIS_T_INVERSE_POWERFUL:	// (p, t^-1 p)
			RingElement::multiplyViaEmbedding(x, a, b);
			x.setBasis(Basis::BASIS_T_INVERSE_POWERFUL);
			break;

		case Basis::BASIS_DECODING:				// (p, d)
			RingElement::multiplyViaEmbedding(x, a, b);
			x.setBasis(Basis::BASIS_DECODING);
			break;

		default:
			return false;
		}
		break;

	case Basis::BASIS_T_INVERSE_POWERFUL:
		switch (b2)
		{
		case Basis::BASIS_T_INVERSE_POWERFUL:	// (t^-1 p, t^-1 p) 
			RingElement::multiplyViaEmbedding(x, a, b);
			x.setBasis(Basis::BASIS_T_INVERSE_POWERFUL);
			break;

		case Basis::BASIS_POWERFUL:				// (t^-1 p, p)
			RingElement::multiplyViaEmbedding(x, a, b);
			x.setBasis(Basis::BASIS_T_INVERSE_POWERFUL);
			break;

		case Basis::BASIS_DECODING:				// (t^-1 p, d)
			b.swapBasesTInversePowerfulAndDecoding(false);
			RingElement::multiplyViaEmbedding(x, a, b);
			x.setBasis(Basis::BASIS_T_INVERSE_POWERFUL);
			break;

		default:
			return false;
		}
		break;

	case Basis::BASIS_DECODING:
		switch (b2)
		{
		case Basis::BASIS_POWERFUL:				// (d, p)
			RingElement::multiplyViaEmbedding(x, a, b);
			x.setBasis(Basis::BASIS_DECODING);
			break;

		case Basis::BASIS_T_INVERSE_POWERFUL:	// (d, t^-1 p)	
			a.swapBasesTInversePowerfulAndDecoding(false);
			RingElement::multiplyViaEmbedding(x, a, b);
			x.setBasis(Basis::BASIS_T_INVERSE_POWERFUL);
			break;
			
		case Basis::BASIS_DECODING:				// (d, d)
			a.swapBasesTInversePowerfulAndDecoding(false);
			b.swapBasesTInversePowerfulAndDecoding(false);
			RingElement::multiplyViaEmbedding(x, a, b);
			x.setBasis(Basis::BASIS_T_INVERSE_POWERFUL);
			break;

		default:
			return false;
		}
		break;
		
	default:
		return false;
	}

	x.setIdealPower(k);
	return true;
}

// Multiplies the two vectors a, b component-wise modulo @param modulus and stores the result in x.
void RingElement::multComponentWiseModulo(coordinate_vec& x, coordinate_vec const& a, coordinate_vec const& b, int modulus)
{
	std::transform(a.cbegin(), a.cend(), b.cbegin(), x.begin(), 
		[modulus](const int& x, const int& y){ return (x * y) % modulus; });
}

// Performs a multiplication via the embedding sigma. That is x = a * b = sigma^-1 ( sigma(a) * sigma(b) ).
// Expects that the bases of a and b are only amongst BASIS_POWERFUL, BASIS_T_INVERSE_POWERFUL and BASIS_DECODING. 
void RingElement::multiplyViaEmbedding(RingElement& x, RingElement const& a, RingElement const& b)
{
	RingLweCryptographyField const& field = a.getField();
	int n = field.getDimension();

	comp_vec sigmaA(n);
	comp_vec sigmaB(n);

	Basis b1 = a.getBasis();
	Basis b2 = b.getBasis();

	// Compute the right transformation matrices to apply sigma in the bases p or d.

	// Default case (p, p)
	TransMat mat1 = TransMat::COMPLEX_CRT_M;
	TransMat mat2 = TransMat::COMPLEX_CRT_M;
	TransMat mat3 = TransMat::COMPLEX_CRT_M_INVERSE;

	if (b1 == Basis::BASIS_POWERFUL || b1 == Basis::BASIS_T_INVERSE_POWERFUL){
		if (b2 == Basis::BASIS_DECODING){			//(p, d)
			mat2 = TransMat::COMPLEX_CRT_M_STAR_INVERSE;
			mat3 = TransMat::COMPLEX_CRT_M_STAR;
		}
	}
	if (b1 == Basis::BASIS_DECODING){
		if (b2 == Basis::BASIS_POWERFUL || b2 == Basis::BASIS_T_INVERSE_POWERFUL){	//(d, p), (d, t^-1 p)
			mat1 = TransMat::COMPLEX_CRT_M_STAR_INVERSE;
			mat3 = TransMat::COMPLEX_CRT_M_STAR;
		}
		if (b2 == Basis::BASIS_DECODING){			//(d, d)
			mat1 = TransMat::COMPLEX_CRT_M_STAR_INVERSE;
			mat2 = TransMat::COMPLEX_CRT_M_STAR_INVERSE;
			mat3 = TransMat::COMPLEX_CRT_M_STAR;
		}
	}

	// convert a and b into complex vecs
	comp_vec aComp = comp_vec(n);
	convertIntVecToCompVec(aComp, a.getCoordinates());
	comp_vec bComp = comp_vec(n);
	convertIntVecToCompVec(bComp, b.getCoordinates());
		
	// embed a and b into H
	field.applyCompTransformation(mat1, sigmaA, aComp);
	field.applyCompTransformation(mat2, sigmaB, bComp);

	// component-wise multiplication in H
	RingElement::multComponentWise<complex_type>(sigmaA, sigmaA, sigmaB);

	// recover x = a * b via the inverse embedding
	field.applyCompTransformation(mat3, aComp, sigmaA);
	convertCompVecToIntVec(x.getCoordinates(), aComp);
}

// Multiplies two elements in the basis CRT. Thus it's just component-wise multiplication.
void RingElement::multiplyInCrt(RingElement& x, RingElement const& a, RingElement const& b)
{
	multComponentWiseModulo(x.getCoordinates(), a.getCoordinates(), b.getCoordinates(), x.getPrimeModulusQ());
}


/******************************************\
		  Special multiplications
\******************************************/

// x = t * a
bool multWithT(RingElement& x, RingElement const& a)
{
	if (x != a){
		x = a;
	}
	Basis b = a.getBasis();
	switch (b)
	{
	case Basis::BASIS_DECODING:
		x.changeBasisTo(Basis::BASIS_T_INVERSE_POWERFUL);
	case Basis::BASIS_T_INVERSE_POWERFUL:
		// t * a = t <t^-1 p, a'> = <p, a'>
		// remove the factor t^-1 from the basis.
		x.setBasis(Basis::BASIS_POWERFUL); 
		x.setIdealPower(0);
		return true;

	case Basis::BASIS_T_INVERSE_CRT:
		// t * a = t <t^-1 c, a'> = <c, a'>
		// remove the factor t^-1 from the basis.
		x.setBasis(Basis::BASIS_CRT);
		x.setIdealPower(0);
		return true;

	default:
		return false;
	}
}

// x = t^-1 * a
bool multWithT_Inverse(RingElement& x, RingElement const& a)
{
	Basis b = a.getBasis();
	x = a;
	
	switch (b)
	{
	case Basis::BASIS_POWERFUL:
		// t^-1 * a = t^-1 <p, a'> = <t^-1 p, a'>
		// add the factor t^-1 to the basis.
		x.setBasis(Basis::BASIS_T_INVERSE_POWERFUL);
		x.changeBasisTo(Basis::BASIS_DECODING);
		x.setIdealPower(1);
		return true;

	case Basis::BASIS_CRT:
		// t^-1 * a = t^-1 <c, a'> = <t^-1 c, a'>
		// add the factor t^-1 to the basis.
		x.setBasis(Basis::BASIS_T_INVERSE_CRT);
		x.setIdealPower(1);
		return true;

	case Basis::BASIS_DECODING:
	case Basis::BASIS_T_INVERSE_POWERFUL:
	case Basis::BASIS_T_INVERSE_CRT:
	{// usual multiplication with t^-1
		x *= a.getField().getElementT_inverse();
	}
		return true;

	default:
		return false;
	}
}

/******************************************\
			Modulo operations
\******************************************/

// x = a mod q. That is, the coordinate vector of a is taken component-wise mod q.
// For an integer e we want 0 <= e mod q< q.
void mod(RingElement& x, RingElement const& a, int q)
{
	mod(x.getCoordinates(), a.getCoordinates(), q);
	x.setBasis(a.getBasis());
	x.setPrimeModulusQ(q);
	x.setIdealPower(a.getIdealPower());
}

// x = a mod q. Take a component-wise mod q. For an integer e we want 0 <= e mod q< q.
void mod(coordinate_vec& x, coordinate_vec const& a, int q)
{
	if (x.size() != a.size()){
		x.resize(a.size());
	}
	std::transform(a.cbegin(), a.cend(), x.begin(), [q](coordinate_type const& coord){ 
		return (coord % q);
	});
	// in c++, the modulo value of a negative number is also negative.
	// we want 0 <= e mod q < q
	std::for_each(x.begin(), x.end(), [q](coordinate_type& e){
		if (e < 0){
			e += q;
		}
	});
}

RingElement operator%(RingElement const& a, int q)
{
	RingElement x(a.getFieldSptr());
	mod(x, a, q);
	return x;
}

coordinate_vec operator%(coordinate_vec const& a, int q)
{
	coordinate_vec x = coordinate_vec();
	mod(x, a, q);
	return x;
}

/************************************************************************/
/*                      decoding                                        */
/************************************************************************/

// decoding procedure from the companion work [May15]
void decode(RingElement& x, RingElement const& a)
{
	if (x != a){
		x = a;
	}
	LweField const& field = x.getField();
	pos_int k = x.getIdealPower();
	coordinate_vec& coords = x.getCoordinates();
	if (k == 0){
		throw "Can't decode elements from R or R_q!";
	}
	if (k == 1){
		if (x.changeBasisTo(Basis::BASIS_DECODING)){
			computeUniqueRepresentative(x, x);
		}
	}
	else{
		switch (x.getBasis())
		{
		case Basis::BASIS_T_INVERSE_CRT:
		{
			RingElement g = field.getElementG();
			RingElement g_pow = g;
			// compute g^(k-1)
			for (pos_int i = 0; i < k - 2; i++)
			{
				g_pow *= g;
			}
			// multiply x with g^(k-1) in the CRT basis
			x *= g_pow;
		}
			break;

		case Basis::BASIS_DECODING:
			// multiply x with g^(k-1) in the decoding basis
			for (pos_int i = 0; i < k - 1; i++)
			{
				field.applyIntTransformation(TransMat::INTEGER_MULT_G_DEC, coords, coords);
				if (x.getPrimeModulusQ() > 0){
					x.setCoordinates(coords % x.getPrimeModulusQ());
				}
			}
			break;

		case Basis::BASIS_T_INVERSE_POWERFUL:
			// multiply x with g^(k-1) in the powerful basis
			for (pos_int i = 0; i < k - 1; i++)
			{
				field.applyIntTransformation(TransMat::INTEGER_MULT_G_POW, coords, coords);
				if (x.getPrimeModulusQ() > 0){
					x.setCoordinates(coords % x.getPrimeModulusQ());
				}
			}
			break;
		}

		if (x.changeBasisTo(Basis::BASIS_DECODING)){
			computeUniqueRepresentative(x, x);
		}

		for (pos_int i = 0; i < k - 1; i++)
		{
			try{// try to divide x by g^(k-1) in the decoding basis.
				field.applyIntTransformation(TransMat::INTEGER_DIV_G_DEC, coords, coords);
			}
			catch (std::exception& e){// decoding failure
				std::cout << "An exception occurred: " << e.what();
				return;
			}
		}
	}
}

// computes the unique representative x in R or R^dual for a in R_q or R^dual_q, resp.
bool computeUniqueRepresentative(RingElement& x, RingElement const& a)
{
	coordinate_vec coords = a.getCoordinates();
	int q = a.getPrimeModulusQ();

	// Check if element is a coset representative
	if (q == 0){
		return false;
	}
	// Elements must not be greater or equal to q
	if (std::any_of(coords.cbegin(), coords.cend(), [&q](coordinate_type const& coord){ return coord >= q; })){
		return false;
	}

	//translate elements from [0, q-1] to [-(q-1)/2, (q-1)/2]
	int i = (q - 1) / 2;
	std::for_each(coords.begin(), coords.end(), [&q, &i](coordinate_type& coord){
		if (coord > i){
			coord -= q;
		}
	});

	x.setCoordinates(coords);
	x.setBasis(a.getBasis());
	x.setPrimeModulusQ(0);

	return true;
}


/******************************************\
		Sampling and Discretization
\******************************************/

// Discretizes y to c + p*L(B). c = @param cosetRepresentative, B = c.basis_, y = <@param realCoordinateVector, B>
// p = @param scalingFactor
void discretize(RingElement& x, real_vec& realCoordinateVector, RingElement const& cosetRepresentative, int const scalingFactor)
{
	double sf = scalingFactor;

	std::random_device rDev;
	std::default_random_engine engine(rDev());

	std::transform(realCoordinateVector.cbegin(), realCoordinateVector.cend(), cosetRepresentative.getCoordinates().cbegin(),
		realCoordinateVector.begin(),
		[sf, &engine](double const& y, int const& c){ // lambda function
		double intpart;

		// cut off integer part
		double p = std::modf((c - y) / sf, &intpart);

		// and make sure that p is in [0,1)
		if (p < 0){
			p = p + 1.0;
		}
		std::bernoulli_distribution distribution(1.0 - p);
		if (distribution(engine)){
			return round(y + sf * p); // return this with probability  1 - p
		}
		else{
			return round(y + sf * (p - 1.0)); // otherwise return this
		}
	});
	
	coordinate_vec intCoords(realCoordinateVector.begin(), realCoordinateVector.end());
	int q = cosetRepresentative.getPrimeModulusQ();
	if (q > 0){
		x.setPrimeModulusQ(q);
		intCoords = intCoords % q;
	}
	x.setCoordinates(intCoords);
	x.setBasis(cosetRepresentative.getBasis());
	x.setIdealPower(cosetRepresentative.getIdealPower());
}

// Sample discretized Gaussians over c + pR^dual, where c = @param cosetRepresentative, c.basis_ = decoding basis and p = @param scalingFactor
void sampleDiscretizedGaussian(RingElement& x, real_type mean, real_type stddev, RingElement const& cosetRepresentative, 
	int const scalingFactor /*= 1*/)
{
	real_vec sample = real_vec(x.getCoordinates().size());

	// sample continuous Gaussian in K_RR in the decoding basis.
	sampleGaussiansInK_RR(sample, x.getField(), mean, stddev);

	// scale if necessary
	if (scalingFactor > 1){
		std::for_each(sample.begin(), sample.end(), [&scalingFactor](real_type& y){ y *= scalingFactor; });
	}

	// discretize to c + pR^dual, where c = cosetRepresentative and p = scalingFactor
	discretize(x, sample, cosetRepresentative, scalingFactor);
}

/****************************************\
	Help functions for sampling
\****************************************/

// Sample coordinate vector for continuous Gaussian over K_RR in the decoding basis. K is given by @param field.
void sampleGaussiansInK_RR(real_vec& x, LweField const& field, real_type mean, real_type stddev)
{
	if (x.size() != field.getDimension())
	{
		x.resize(field.getDimension());
	}

	std::random_device rDev;
	std::default_random_engine engine(rDev());
	double m = field.getM(), rad_m = field.getRadM();
	double s = stddev * (sqrt(m / rad_m));
	std::normal_distribution<real_type> gaussian = std::normal_distribution<real_type>(mean, s);

	// sample continuous Gaussian in RR
	std::generate(x.begin(), x.end(), [&engine, &gaussian](){ return gaussian(engine); });

	// transform into Gaussian in H, i.e., x is Gaussian distributed coordinate vector for an Element in K_RR
	field.applyRealTransformation(TransMat::REAL_SAMPLE_GAUSS_D, x, x);
}

// Sample discrete Gaussian over ZZ centered at @param center and with standard deviation @param stddev.
coordinate_type sampleDiscreteGaussianInZZ(real_type stddev, real_type center, real_type funcValue /*= 3*/)
{
	using namespace std;

	random_device rDev;
	default_random_engine engine(rDev());
	
	// compute lower and upper bound
	coordinate_type lb, ub;
	lb = ceil(center - stddev*funcValue);
	ub = floor(center + stddev*funcValue);

	// sample until not rejected
	while (true){ 

		// choose integer uniform at random from [lb, ub] 
		uniform_int_distribution<coordinate_type> uniDist = uniform_int_distribution<coordinate_type>(lb, ub);
		coordinate_type sampleZZ = uniDist(engine);

		// Probability function of Gaussian distribution
		real_type p = exp(-(M_PI * pow(abs(sampleZZ - center), 2) / pow(stddev, 2)));
		bernoulli_distribution bernDist = bernoulli_distribution(p);

		// return with probability p, otherwise reject and repeat
		if (bernDist(engine)){
			return sampleZZ;
		}
	}
}

// Sample discrete Gaussians of standard deviation @param stddev over c+R.
void sampleDiscreteGaussianInR(RingElement& x, RingElement const& c, real_type stddev)
{
	LweField const& field = x.getField();
	int n = field.getDimension();
	real_type log_n = std::log(n);
	real_type center, s;

	coordinate_vec coords = c.getCoordinates();
	real_vec unit_vec = real_vec(n, REAL_TYPE_ZERO);
	real_vec i_th_row_U = real_vec(n);
	real_vec i_th_column_D = real_vec(n);

	// init n-th unit vector 
	unit_vec[n - 1] = REAL_TYPE_ONE;

	for (int i = n-1; i >= 0; i--)
	{
		// compute i-th row of U
		field.applyRealTransformation(TransMat::REAL_GS_DECOMP_U, i_th_row_U, unit_vec);

		// inner product of coords and i-th row of U
		center = std::inner_product(coords.begin(), coords.end(), i_th_row_U.begin(), REAL_TYPE_ZERO);

		field.applyRealTransformation(TransMat::REAL_GS_DECOMP_D, i_th_column_D, unit_vec);
		s = stddev / i_th_column_D[i]; // i-th diagonal entry of D

		// sample discrete Gaussian in the computed range and update coords
		coordinate_type z = sampleDiscreteGaussianInZZ(s, center, log_n);
		coords[i] -= z;

		// init next unit vector
		if (i != 0){
			unit_vec[i] = REAL_TYPE_ZERO;
			unit_vec[i - 1] = REAL_TYPE_ONE;
		}
	}

	x.setCoordinates(coords);
	x.setBasis(Basis::BASIS_POWERFUL);
	x.setPrimeModulusQ(0);
	x.setIdealPower(0);
}

// Samples uniformly distributed entries over [lb,ub] for the coordinate vector of x.
void sampleUniformlyCoordVec(RingElement& x, coordinate_type lb, coordinate_type ub)
{
	using namespace std;

	random_device rDev;
	default_random_engine engine(rDev());
	uniform_int_distribution<coordinate_type> distribution(lb, ub);

	coordinate_vec& coords = x.getCoordinates();
	generate(coords.begin(), coords.end(), [&distribution, &engine](){ return distribution(engine); });
}

/***********************************************\
    Help functions for basis transformations
\***********************************************/

// Change the basis to @param newBasis if possible. Return true iff conversion succeeded.
bool RingElement::changeBasisTo(Basis newBasis){
	if (getBasis() == newBasis){
		return true;
	}
	switch (getBasis())
	{
	case Basis::BASIS_POWERFUL:
		switch (newBasis)
		{
		case Basis::BASIS_T_INVERSE_POWERFUL:
			throw "Can't change basis from p to t^-1p!";
			return false;
		case Basis::BASIS_CRT:
			return swapBasesPowerfulAndCrt(true);
		case Basis::BASIS_T_INVERSE_CRT:
			throw "Can't change basis from p to t^-1c!";
			return false;
		case Basis::BASIS_DECODING:
			throw "Can't change basis from p to d!";
			return false;
		}
		break;

	case Basis::BASIS_T_INVERSE_POWERFUL:
		switch (newBasis)
		{
		case Basis::BASIS_POWERFUL:
			throw "Can't change basis from t^-1p to p!";
			return false;
		case Basis::BASIS_CRT:
			throw "Can't change basis from t^-1p to c!";
			return false;
		case Basis::BASIS_T_INVERSE_CRT:
			return swapBasesPowerfulAndCrt(true);
		case Basis::BASIS_DECODING:
			return swapBasesTInversePowerfulAndDecoding(true);
		}
		break;
		
	case Basis::BASIS_CRT:
		switch (newBasis)
		{
		case Basis::BASIS_POWERFUL:
			return swapBasesPowerfulAndCrt(false);
		case Basis::BASIS_T_INVERSE_POWERFUL:
			throw "Can't change basis from c to t^-1p!";
			return false;
		case Basis::BASIS_T_INVERSE_CRT:
			throw "Can't change basis from c to t^-1c!";
			return false;
		case Basis::BASIS_DECODING:
			throw "Can't change basis from c to d!";
			return false;
		}
		break;
		
	case Basis::BASIS_T_INVERSE_CRT:
		switch (newBasis)
		{
		case Basis::BASIS_POWERFUL:
			throw "Can't change basis from t^-1c to p!";
			return false;
		case Basis::BASIS_T_INVERSE_POWERFUL:
			return swapBasesPowerfulAndCrt(false);
		case Basis::BASIS_CRT:
			throw "Can't change basis from t^-1c to c!";
			return false;
		case Basis::BASIS_DECODING:
			return swapBasesPowerfulAndCrt(false) && swapBasesTInversePowerfulAndDecoding(true);
		}
		break;

	case Basis::BASIS_DECODING:
		switch (newBasis)
		{
		case Basis::BASIS_POWERFUL:
			throw "Can't change basis from d to p!";
			return false;
		case Basis::BASIS_T_INVERSE_POWERFUL:
			return swapBasesTInversePowerfulAndDecoding(false);
		case Basis::BASIS_CRT:
			throw "Can't change basis from d to c!";
			return false;
		case Basis::BASIS_T_INVERSE_CRT:
			return swapBasesTInversePowerfulAndDecoding(false) && swapBasesPowerfulAndCrt(true);
		}
		break;
	}
}

// Converts between the powerful and CRT basis of R and R^dual
bool RingElement::swapBasesPowerfulAndCrt(bool direction)
{
	LweField const& field = getField();
	coordinate_vec& vec = getCoordinates();
	Basis b = getBasis();
	if (direction){ // powerful -> CRT
		if (b != Basis::BASIS_POWERFUL && b != Basis::BASIS_T_INVERSE_POWERFUL){
			return false;
		}
		field.applyIntTransformation(TransMat::INTEGER_CRT_MQ, vec, vec);
		if (b == Basis::BASIS_POWERFUL){ // p -> c
			setBasis(Basis::BASIS_CRT);
		}
		else{ // t^-1 p -> t^-1 c
			setBasis(Basis::BASIS_T_INVERSE_CRT);
		}
		// modulus of element might still be zero
		setPrimeModulusQ(field.getModulus());
		return true;
	}
	else{ // CRT -> powerful
		if (b != Basis::BASIS_CRT && b != Basis::BASIS_T_INVERSE_CRT){
			return false;
		}
		field.applyIntTransformation(TransMat::INTEGER_CRT_MQ_INVERSE, vec, vec);
		if (b == Basis::BASIS_CRT){ // c -> p
			setBasis(Basis::BASIS_POWERFUL);
		}
		else{ // t^-1 c -> t^-1 p
			setBasis(Basis::BASIS_T_INVERSE_POWERFUL);
		}
		
		return true;
	}
}

// Converts between the powerful and decoding basis of R^dual
bool RingElement::swapBasesTInversePowerfulAndDecoding(bool direction)
{
	LweField const& field = getField();
	coordinate_vec& vec = getCoordinates();
	Basis b = getBasis();
	if (direction){ // powerful -> decoding
		if (b != Basis::BASIS_T_INVERSE_POWERFUL){
			return false;
		}
		field.applyIntTransformation(TransMat::INTEGER_L_M_INVERSE, vec, vec);
		setBasis(Basis::BASIS_DECODING);
	}
	else{ // decoding -> powerful
		if (b != Basis::BASIS_DECODING){
			return false;
		}
		field.applyIntTransformation(TransMat::INTEGER_L_M, vec, vec);
		setBasis(Basis::BASIS_T_INVERSE_POWERFUL);
	}
	int q = getPrimeModulusQ();
	if (q > 0){ // compute manually modulo q, since L_m is a normal matrix over ZZ
		vec = vec % q;
	}
	return true;
}


void convertIntVecToCompVec(comp_vec& compVec, coordinate_vec const& intVec)
{
	if (compVec.size() < intVec.size()){
		compVec.resize(intVec.size());
	}
	std::transform(intVec.cbegin(), intVec.cend(), compVec.begin(), [](int value){ return complex_type(value); });
}

void convertCompVecToIntVec(coordinate_vec& intVec, comp_vec const& compVec)
{
	if (intVec.size() < compVec.size()){
		intVec.resize(compVec.size());
	}
	std::transform(compVec.cbegin(), compVec.cend(), intVec.begin(), [](complex_type value){ return std::round(value.real()); });
}

/********************\
		Output
\********************/

std::ostream& operator<<(std::ostream& out, RingLweCryptographyElement const& elmt)
{
	using namespace std;
	string basis;
	switch (elmt.getBasis())
	{
	case Basis::BASIS_POWERFUL:
		basis = "POWERFUL";
		break;
	case Basis::BASIS_CRT:
		basis = "CRT";
		break;
	case Basis::BASIS_DECODING:
		basis = "DECODING";
		break;
	case Basis::BASIS_T_INVERSE_POWERFUL:
		basis = "t^-1 POWERFUL";
		break;
	case Basis::BASIS_T_INVERSE_CRT:
		basis = "t^-1 CRT";
		break;
	}
	int n = elmt.getField().getDimension();
	out << "Dimension:" << endl << n << endl << endl;
	out << "Basis:" << endl << basis << endl << endl;
	pos_int k = elmt.getIdealPower();
	if (k > 1){
		out << "Ideal power: " << endl << k << endl << endl;
	}
	out << "Coordinates:" << endl;
	coordinate_vec vec = elmt.getCoordinates();
	out << "(";
	for (int i = 0; i < n; i++)
	{
		if (i == n - 1){
			out << vec[i];
		}
		else{
			out << vec[i] << ",";
		}
	}
	out << ")" << endl << endl;

	pos_int q = elmt.getPrimeModulusQ();
	if (q > 0){
		out << "Modulus:" << endl << q << endl << endl;
	}
	return out;
}







