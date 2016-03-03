#pragma once

#ifndef DSC_H
#define DSC_H

#include <stdio.h>
#include "Example util.h"
#include "Main_Classes/RingLweCryptographyField.h"

typedef RingElement::Basis Basis;
typedef std::vector<RingElement> elmt_vec;

using namespace RLWE_Toolkit;

/*
Example implementation of the "Dual-Style Cryptosystem" from the companion work [May15].
Shows how the toolkit can be used.
*/
class DualStyleCryptosystem
{
public:
	DualStyleCryptosystem(int m, int p, int q, int l);

	~DualStyleCryptosystem();

	// Key generation
	void gen();

	// Encryption
	elmt_vec enc(RingElement& message);

	// Decryption
	RingElement dec(elmt_vec& ciphertext);

	elmt_vec public_key_;

	// The cyclotomic field K
	std::shared_ptr<RingLweCryptographyField> field_;

private:

	elmt_vec secret_key_;

	int p_; // Represents the message space R_p
	int l_; // Keys and cyphertexts are vectors of length l

	// standard deviation for the discrete Gaussian distribution over R
	real_type discGaussian_;

	// standard deviation for the Gaussian LWE error distribution
	real_type LweErrorDist_; 
};

DualStyleCryptosystem::DualStyleCryptosystem(int m, int p, int q, int l) :
p_(p),
l_(l)
{
	field_ = std::make_shared<RingLweCryptographyField>(m, q);
	int n = field_->getDimension();

	discGaussian_ = sqrt(n) * log(n); 
	LweErrorDist_ = (log(n) / q);	

	gen();
}

DualStyleCryptosystem::~DualStyleCryptosystem(){}

void DualStyleCryptosystem::gen()
{
	RingElement ph = RingElement(field_), zero = field_->getZero(RingElement::Basis::BASIS_POWERFUL);
	elmt_vec a = elmt_vec(l_ + 1, ph);
	elmt_vec x = elmt_vec(l_ + 1, ph);
	int q = field_->getModulus();

	// a_0 = -1
	RingElement a0 = -field_->getOne(Basis::BASIS_CRT);
	a[0] = a0;

	// x_0 <- D_r
	RingElement x0 = RingElement(field_);
	sampleDiscreteGaussianInR(x0, zero, discGaussian_);
	x[0] = x0;

	for (int i = 1; i < l_; i++){

		// sample a_1 to a_(l-1) independent and uniformaly at random
		sampleUniformlyCoordVec(a0, 0, q - 1);
		a[i] = a0;

		// sample x_1 to x_(l-1) from the discrete Gaussian distribution D_r
		sampleDiscreteGaussianInR(x0, zero, discGaussian_);
		x[i] = x0;
	}

	// a_l = - <a', x'> where a' = (a_0, ..., a_(l-1)) and x' = (x_0, ..., x_(l-1))
	a0.setCoordinates(coordinate_vec(field_->getDimension(), 0));
	for (int i = 0; i < l_; i++){
		a0 += (a[i] * x[i]);
	}
	a[l_] = -a0;

	// x_l = 1
	x[l_] = field_->getOne(Basis::BASIS_POWERFUL);

	// a = (a_1, ..., a_l) is public key,
	public_key_ = elmt_vec(a.begin() + 1, a.end());

	// x = (x_1, ..., x_l) is secret key;
	secret_key_ = elmt_vec(x.begin() + 1, x.end());
}

elmt_vec DualStyleCryptosystem::enc(RingElement& message)
{
	RingElement ph = RingElement(field_);
	elmt_vec errors = elmt_vec(l_ + 1, ph);
	elmt_vec cipher = elmt_vec(l_, ph);

	// sample short errors
	RingElement e0 = RingElement(field_);
	RingElement zero = field_->getZero(Basis::BASIS_DECODING);
	std::generate(errors.begin(), errors.end() - 1, [&](){
		sampleDiscretizedGaussian(e0, REAL_TYPE_ZERO, sqrt(LweErrorDist_), zero, p_);
		return e0;
	});

	// hide message in e_l
	// discretize sample to t^-1 * message + pR^dual
	RingElement cosetRep = RingElement(field_);
	multWithT_Inverse(cosetRep, message);
	sampleDiscretizedGaussian(e0, REAL_TYPE_ZERO, sqrt(LweErrorDist_), cosetRep, p_);
	errors[l_] = e0;

	// encrypt with public key
	for (int i = 0; i < l_; i++){
		cipher[i] = (errors[0] * public_key_[i]) + errors[i + 1];
	}

	return cipher;
}

RingElement DualStyleCryptosystem::dec(elmt_vec& ciphertext)
{
	// Decrypt with secret key
	// decode d = <c,x> (c = ciphertext, x = secret key)
	RingElement d = field_->getZero(Basis::BASIS_DECODING);
	for (int i = 0; i < l_; i++){
		d += (ciphertext[i] * secret_key_[i]);
	}
	decode(d,d);

	// message = t*d mod pR
	RingElement message = RingElement(field_);
	multWithT(message, d);
	message = message % p_;

	return message;
}

#endif // !DSC_H