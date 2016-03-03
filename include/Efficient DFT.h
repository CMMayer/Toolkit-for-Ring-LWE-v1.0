#pragma once

#ifndef EDFT_H

#define EDFT_H

#include <complex>
#include <vector>

#include "Data_Type_Settings.h"

typedef std::complex<real_type> complex_type;
typedef std::vector<complex_type> complex_vec;

namespace RLWE_Toolkit {
	namespace DFT {

		/*
		Cooley-Tukey FFT for power of two input.size() = N = 2^k.
		Calls cooley_tukey_fft(output, input, N)
		*/
		void cooley_tukey_fft(complex_vec& output, complex_vec const& input);

		/*
		Cooley-Tukey IFFT for power of two input.size() = N = 2^k.
		Calls cooley_tukey_ifft(output, input, N)
		*/
		void cooley_tukey_ifft(complex_vec& output, complex_vec const& input);

		/*
		Recursive implementation of Cooley-Tukey FFT for powers of two N = 2^k.
		The whole recursion works on the input vector. The values @param start and @param step
		are used to traverse through certain sub vectors throughout the recursion.
		*/
		void cooley_tukey_fft(complex_vec& output, complex_vec const& input, int N, int start = 0, int step = 1);

		/*
		Recursive implementation of Cooley-Tukey IFFT for powers of two N = 2^k.
		The whole recursion works on the input vector. The values @param start and @param step
		are used to traverse through certain sub vectors throughout the recursion.
		*/
		void cooley_tukey_ifft(complex_vec& output, complex_vec const& input, int N, int start = 0, int step = 1);

		/*
		Implementation of Rader's DFT for prime length input. For a detailed explanation see:
		http://www.skybluetrades.net/blog/posts/2013/12/31/data-analysis-fft-9.html
		*/
		void rader_dft_for_primes(complex_vec& output, complex_vec const& input, complex_vec const& precomp_DFT_omega_p, int generator);

		/*
		Slightly adjusted implementation of Rader's DFT to fit an application of the Chinese remainder transform.
		*/
		void rader_crt_for_primes(complex_vec& output, complex_vec const& input, complex_vec const& precomp_DFT_omega_p, int generator);

		/*
		Slightly adjusted implementation of Rader's DFT to fit an application of the adjoint Chinese remainder transform.
		*/
		void rader_crt_star_for_primes(complex_vec& output, complex_vec const& input, complex_vec const& precomp_DFT_omega_p, int generator);

		/*
		Applies the permutation induced by the generator @param generator of the unit group of ZZ_p.
		Input might have size p or (p-1), output must have size (p-1).
		*/
		void unitsGeneratorPermutation(complex_vec& output, complex_vec const& input, int generator, int p);

		/*
		Applies the permutation induced by the inverse of the generator @param generator of the unit group of ZZ_p.
		Input must have size (p-1), Output might have size p or (p-1).
		*/
		void unitsGeneratorPermutationInverse(complex_vec& output, complex_vec const& input, int generator, int p);

		/*
		Zero pads the input vector to the new size. @param newSize must be >= 2*input.size() - 3.
		zeroPadding inserts the needed amount of zeros between the first and second entry of the input vector.
		*/
		void zeroPadding(complex_vec& output, complex_vec const& input, int newSize);

		/*
		Creates the output vector of size @param newSize, whose entries are a cyclic padding of the input.
		For example, (a b c) -> (a b c a b c a b). @param newSize must be >= 2*input.size() - 3.
		*/
		void cyclicPadding(complex_vec& output, complex_vec const& input, int newSize);

	}
}

#endif // !EDFT_H


