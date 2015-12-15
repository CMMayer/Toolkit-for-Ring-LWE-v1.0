#pragma once

#ifndef CPKC_H
#define CPKC_H

#include <stdio.h>
#include "Example util.h"
#include "RLWE/RingLweCryptographyField.h"

#define REAL_TYPE_ZERO 0.0
#define REAL_TYPE_ONE 1.0

typedef std::pair<RingElement*, RingElement*> elmt_pair;
typedef RingElement::Basis Basis;
typedef unsigned int pos_int;

/*
Example implementation of the "Compact Public-Key Cryptosystem" from the companion work [May15].
Shows how the toolkit can be used.
*/
class CompactPublicKeyCryptosystem
{
public:
	CompactPublicKeyCryptosystem(pos_int m, pos_int p, pos_int q);

	~CompactPublicKeyCryptosystem();

	// Key generation
	void gen();

	// Encryption
	elmt_pair enc(RingElement& message);

	// Decryption
	RingElement dec(elmt_pair& ciphertext);

	elmt_pair public_key_;

	// The cyclotomic field K
	std::shared_ptr<RingLweCryptographyField> field_;

private:

	RingElement *secret_key_;

	int p_; // Defines the message space R_p

	// standard deviation for the Gaussian LWE error distribution
	real_type LweErrorDist_; 
};

CompactPublicKeyCryptosystem::CompactPublicKeyCryptosystem(pos_int m, pos_int p, pos_int q) :
p_(p)
{
	field_ = std::make_shared<RingLweCryptographyField>(m, q);
	LweErrorDist_ = log(field_->getDimension()) / q;
	gen();
}

CompactPublicKeyCryptosystem::~CompactPublicKeyCryptosystem()
{
	delete secret_key_;
	secret_key_ = 0;
}

void CompactPublicKeyCryptosystem::gen()
{
	int n = field_->getDimension();
	int q = field_->getModulus();
	RingElement *a = new RingElement(field_, Basis::BASIS_CRT, coordinate_vec(n, 0), q), e = RingElement(field_), zero = field_->getZero(Basis::BASIS_DECODING);
	secret_key_ = new RingElement(field_);

	// sample uniformly a, and x (=secret key) and e from the discretized error dist
	sampleUniformlyCoordVec(*a, 0, q - 1);
	sampleDiscretizedGaussian(*secret_key_, REAL_TYPE_ZERO, sqrt(LweErrorDist_), zero);
	sampleDiscretizedGaussian(e, REAL_TYPE_ZERO, sqrt(LweErrorDist_), zero, p_);

	// b = t*g*(a*x + e)
	RingElement *b = new RingElement(field_->getElementG()*((*a) * (*secret_key_) + e));
	multWithT(*b, *b);

	public_key_ = elmt_pair(a, b);
}

elmt_pair CompactPublicKeyCryptosystem::enc(RingElement& message)
{
	RingElement z = RingElement(field_), e1 = RingElement(field_), e2 = RingElement(field_), cosetRep = RingElement(field_), zero = field_->getZero(Basis::BASIS_DECODING);

	// sample z and e'
	sampleDiscretizedGaussian(z, REAL_TYPE_ZERO, sqrt(LweErrorDist_), zero);
	sampleDiscretizedGaussian(e1, REAL_TYPE_ZERO, sqrt(LweErrorDist_), zero, p_);

	// hide message in e''
	// sample e'' discretized to t^-1 * message + pR^dual
	multWithT_Inverse(cosetRep, message);
	sampleDiscretizedGaussian(e2, REAL_TYPE_ZERO, sqrt(LweErrorDist_), cosetRep, p_);

	// Encrypt with public key
	// u = t*g*(a*z + e'), v = z*b + e''
	RingElement *u = new RingElement(field_->getElementG()*(z * (*public_key_.first) + e1));
	multWithT(*u, *u);
	RingElement *v = new RingElement((z * (*public_key_.second) + e2));

	// (u,v) is the ciphertext
	return elmt_pair(u, v);
}

RingElement CompactPublicKeyCryptosystem::dec(elmt_pair& ciphertext)
{
	// Decrypt with secret key
	// decode d = v - u*x
	RingElement d = *ciphertext.second - *ciphertext.first * *secret_key_;
	decode(d,d);

	// message = t*d mod pR
	RingElement message = RingElement(field_);
	multWithT(message, d);
	message = message % p_;

	return message;
}

#endif // !CPKC_H
