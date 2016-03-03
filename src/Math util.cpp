#define _USE_MATH_DEFINES

#include <math.h>
#include "boost/math/special_functions/prime.hpp"
#include "boost/math/common_factor.hpp"
#include "boost/multiprecision/miller_rabin.hpp"

#include <random>
#include <set>

#include "Math util.h"

// Finds a generator of the units of ZZ_p.
int RLWE_Toolkit::Math_util::findGeneratorOfZZpUnits(int p)
{
	using namespace std;

	// compute prime power factorization of phi(p) = p-1
	std::list<std::pair<pos_int, pos_int>> primePowers = primeFactorization(p - 1);

	random_device rDev;
	default_random_engine engine(rDev());
	uniform_int_distribution<int> distribution(1, p-1);
	set<int> test;
	while (true){
		// uniformly random element in [1, p-1]
		int g = distribution(engine);
		if (test.count(g) != 0){ // test if g was already checked
			continue;
		}
		bool flag = true;
		for each (auto& pair in primePowers)
		{
			int q = pair.first;
			// Look if g^((p-1)/q) != 1 for each prime divisor q of p-1.
			if (fastModPow(g, (p - 1) / q, p) == 1){
				flag = false;
				test.insert(g);
				break;
			}
		}
		if (flag){
			return g;
		}
	}
}

// Fast exponentiation modulo @param modulus.
int RLWE_Toolkit::Math_util::fastModPow(int base, int exponent, int modulus)
{
	int result = 1;
	// use the binary representation of the exponent
	while (exponent > 0){
		if (exponent % 2 != 0){
			result = (result*base) % modulus;
		}
		base = (base*base) % modulus;
		exponent /= 2;
	}
	return result;
}

// Computes the prime factorization of m. Output is given by a list of pairs. For each prime p dividing m a pair
// (p, k) contains the prime factor and the prime exponent, i.e., p^k|m but not p^k+1|m.
// This algorithm works only for numbers m whose prime factors are within the first 10000 primes.
std::list<std::pair<pos_int, pos_int>> RLWE_Toolkit::Math_util::primeFactorization(pos_int const m)
{
	using namespace boost::math;
	pos_int n = m;

	std::list<std::pair<pos_int, pos_int>> result = std::list<std::pair<pos_int, pos_int>>();
	for (int i = 0; i < max_prime; i++){ //max_prime = 10000, i.e., the 10000-th prime is the largest prime we can test
		if (n == 1){
			break;
		}
		// i-th prime number (list starts at 2)
		pos_int p = prime(i);
		pos_int k = 0;

		// count how often p divides m
		while ((n % p) == 0){
			k++;
			n /= p;
		}
		if (k > 0){
			std::pair<pos_int, pos_int> pair(p, k);
			result.push_back(pair);
		}
	}
	return result;
}

/*
* Euler's totient function phi(n).
* http://en.wikipedia.org/wiki/Euler%27s_totient_function
*
* This is an *EXTREMELY* fast function and uses
* several tricks to recurse.
*
* It assumes you have a list of primes and a fast
* isprime() function.  Typically, you use a bitset
* to implement the sieve of Eratosthenes and use
* isprime() on that.  You should also have a vector
* of all known primes below a certain limit.
*
* Additionally, you should have a fast gcd algorithm.
*
* So, we have three dependencies here:
*
* - isprime(int) (typically look up bitset sieve)
* - coordinate_vec primes (vector of prime numbers, typically sieved)
* - binary_gcd(int, int) or similar, fast, gcd function.
*
* This function is placed in the public domain by the author,
* Christian Stigen Larsen, http://csl.sublevel3.org
*
*********************************************************************************
*
* The dependencies are implemented via the Boost-library.
*
*/
pos_int const RLWE_Toolkit::Math_util::eulerTotient(pos_int const n)
{
	using namespace boost::math;

	// Base case
	if (n < 2)
		return 0;

	// Lehmer's conjecture
	if ((n == 2) || boost::multiprecision::miller_rabin_test(n, 25))
		return n - 1;

	// Even number?
	if ((n & 1) == 0) {
		pos_int m = n >> 1;
		return !(m & 1) ? eulerTotient(m) << 1 : eulerTotient(m);
	}

	// For all primes ...
	for (int i = 0; i < max_prime; i++) //traverse through primes up to the 10000-th prime
	{
		pos_int m = prime(i);
		if (n % m) continue;

		// phi is multiplicative for coprime integers
		pos_int o = n / m;
		pos_int d = gcd(m, o);
		return d == 1 ? eulerTotient(m)*eulerTotient(o) : eulerTotient(m)*eulerTotient(o)*d / eulerTotient(d);
	}
}

// Searches the first prime p > min, s.t. p = 1 mod m.
pos_int RLWE_Toolkit::Math_util::getPrimeMod1(pos_int const m, pos_int const min)
{
	using namespace boost::math;
	int i = 1;
	pos_int n = prime(i++);
	while (n < min){ // search the first prime greater than min
		n = prime(i++);
	}
	while ((n % m) != 1){ // now search for the first prime congruent 1 mod m
		n = prime(i++);
	}

	return n;
}

int RLWE_Toolkit::Math_util::getNextPowerOfTwo(int n)
{
	if (n <= 0){
		return 0;
	}

	// check if n is already a power of two
	if (isPowerOfTwo(n)){
		return n;
	}

	int res = 2;
	// shift bits of n to the right until n = 0.
	while ((n = n >> 1)){
		res *= 2;
	}
	return res;
}

// Returns true iff n is a power of two
bool RLWE_Toolkit::Math_util::isPowerOfTwo(int n)
{
	// (n-1) & n is zero iff n is power of two
	return !((n - 1) & n);
}

// Computes the value of the "first" m-th root of unity in CC.
complex_type const RLWE_Toolkit::Math_util::computeRootOfUnity(pos_int m)
{
	real_type m_temp = m;
	// COMP_M_I = sqrt(-1)
	return exp(2 * M_PI * COMP_M_I / m_temp);
}

// Searches for an element of order @param order. Succeeds if @param order divides p-1
NTL::ZZ_p const RLWE_Toolkit::Math_util::findElementOfOrder(pos_int order, pos_int p)
{
	using namespace NTL;
	ZZ q = ZZ(p);
	ZZ_p::init(q);

	for (ZZ_p g = ZZ_p(2); rep(g) < q; g++){

		// is g a unit in ZZq, i.e., gcd(g,p) = 1?
		if (IsOne(GCD(rep(g), q))){

			// g^order = 1?
			if (IsOne(power(g, order))){

				bool flag = true;

				// check if g^(order / p) != 1 for all primes p|order
				for each (std::pair<pos_int, pos_int> pair in primeFactorization(order))
				{
					pos_int p = pair.first;
					if (IsOne(power(g, order / p))) {
						flag = false;
						break;
					}
				}
				if (flag){
					return g;
				}
			}
		}
	}
}