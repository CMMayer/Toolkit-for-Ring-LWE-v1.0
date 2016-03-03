#pragma once

#ifndef MU_H
#define MU_H

#include <list>
#include <tuple>
#include <vector>
#include <complex>

#include "Data_Type_Settings.h"

#include "NTL/ZZ_p.h"

namespace RLWE_Toolkit {
	namespace Math_util {

		// For use in FFT algorithms 

		// compute a generator of the units of ZZ_p.
		int findGeneratorOfZZpUnits(int p);

		// compute a power of two 2^k >= n
		int getNextPowerOfTwo(int n);

		// true iff n is a power of two, n = 2^k.
		bool isPowerOfTwo(int n);

		// Fast exponentiation modulo @param modulus. Return base^exponent mod modulus.
		int fastModPow(int base, int exponent, int modulus);

		// Returns phi(n) (Euler totient)
		pos_int const eulerTotient(pos_int const n);

		// Returns a list of pairs. For each prime p with p^k|m and not p^(k+1)|m there is a pair (p, k).
		std::list<std::pair<pos_int, pos_int>> primeFactorization(pos_int const m);

		// Returns the first prime q >= min, s.t. q = 1 mod m.
		pos_int getPrimeMod1(pos_int const m, pos_int const min);

		// Computes the complex value of the "first" m-th root of unity. First in the sense,
		// that roots of unity are ordered counterclockwise in a circle in the complex plane.
		complex_type const computeRootOfUnity(pos_int m);

		/*
		* Searches an element of order @param order in the multiplicative group ZZ_p^*.
		* The order of the group ZZ_p^* is order(ZZ_p^*) = eulerTotient(p) = p-1.
		* For each divisor d of order(ZZ_p^*), ZZ_p^* contains exactly eulerTotient(d)
		* many distinct elements of order d. Hence the procedure will successfully
		* return an element of the desired order, if this order divides p-1.
		*
		* @param order: Order of the element to find.
		*
		* @return Return an element of order @param order, if order divides m.
		*/
		NTL::ZZ_p const findElementOfOrder(pos_int order, pos_int p);

	}
}

#endif // !MU_H