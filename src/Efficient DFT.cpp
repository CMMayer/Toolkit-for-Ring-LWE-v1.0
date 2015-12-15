#define _USE_MATH_DEFINES

#include <math.h>
#include <algorithm>
#include <numeric>

#include "Math util.h"

#include "Efficient DFT.h"

#define REAL_TYPE_ZERO real_type(0)
#define REAL_TYPE_ONE real_type(1)
#define COMP_M_I complex_type(REAL_TYPE_ZERO, REAL_TYPE_ONE)
#define COMP_M_SQRT1_2 complex_type(M_SQRT1_2)

/*
Implementation of Rader's DFT for prime length input. For a detailed explanation see:
http://www.skybluetrades.net/blog/posts/2013/12/31/data-analysis-fft-9.html
*/
void rader_dft_for_primes(complex_vec& output, complex_vec const& input, complex_vec const& precomp_DFT_omega_p, int generator)
{
	int p = input.size();

	complex_vec permuted_vec = complex_vec(p - 1), zero_padded_vec;

	complex_type h0 = input[0];

	// the first entry of the output is just the sum of all input entries.
	output[0] = std::accumulate(input.begin(), input.end(), complex_type(REAL_TYPE_ZERO));

	// apply the permutation given by the generator g.
	unitsGeneratorPermutation(permuted_vec, input, generator, p);

	// zero pad to a vector of size equal to the next power of two if necessary.
	if (!isPowerOfTwo(p - 1)){
		zeroPadding(zero_padded_vec, permuted_vec, precomp_DFT_omega_p.size());
	}
	else{
		zero_padded_vec = permuted_vec;
	}

	// compute DFT of this vector via Cooley-Tukey FFT.
	cooley_tukey_fft(zero_padded_vec, zero_padded_vec);

	// component-wise multiply with the pre-computed DFT of the vector of powers of the p-th root of unity.
	std::transform(zero_padded_vec.begin(), zero_padded_vec.end(), precomp_DFT_omega_p.begin(), zero_padded_vec.begin(), std::multiplies<complex_type>());

	// compute IFFT via Cooley-Tukey.
	cooley_tukey_ifft(zero_padded_vec, zero_padded_vec);

	// extract result, i.e., the first p - 1 entries.
	std::copy(zero_padded_vec.begin(), zero_padded_vec.begin() + (p - 1), permuted_vec.begin());

	// undo the permutation.
	unitsGeneratorPermutationInverse(output, permuted_vec, generator, p);

	// add h0.
	std::for_each(output.begin() + 1, output.end(), [h0](complex_type& entry){ entry += h0; });
}

/*
Slightly changed implementation of Rader's DFT for prime length input. This algorithm computes the CRT instead of the DFT.
*/
void rader_crt_for_primes(complex_vec& output, complex_vec const& input, complex_vec const& precomp_DFT_omega_p, int generator)
{
	int p = input.size() + 1;

	complex_vec adjusted_input = complex_vec(p), permuted_vec = complex_vec(p - 1), zero_padded_vec;

	std::copy(input.begin(), input.end(), adjusted_input.begin());
	adjusted_input[p - 1] = complex_type(REAL_TYPE_ZERO);

	complex_type h0 = adjusted_input[0];

	// apply the permutation given by the generator g.
	unitsGeneratorPermutation(permuted_vec, adjusted_input, generator, p);

	// zero pad to a vector of size equal to the next power of two if necessary.
	if (!isPowerOfTwo(p - 1)){
		zeroPadding(zero_padded_vec, permuted_vec, precomp_DFT_omega_p.size());
	}
	else{
		zero_padded_vec = permuted_vec;
	}

	// compute DFT of this vector via Cooley-Tukey FFT.
	cooley_tukey_fft(zero_padded_vec, zero_padded_vec);

	// component-wise multiply with the pre-computed DFT of the vector of powers of the p-th root of unity.
	std::transform(zero_padded_vec.begin(), zero_padded_vec.end(), precomp_DFT_omega_p.begin(), zero_padded_vec.begin(), std::multiplies<complex_type>());

	// compute IFFT via Cooley-Tukey.
	cooley_tukey_ifft(zero_padded_vec, zero_padded_vec);

	// extract result, i.e., the first p - 1 entries.
	std::copy(zero_padded_vec.begin(), zero_padded_vec.begin() + (p - 1), permuted_vec.begin());

	// undo the permutation.
	unitsGeneratorPermutationInverse(output, permuted_vec, generator, p);

	// add h0.
	std::for_each(output.begin(), output.end(), [h0](complex_type& entry){ entry += h0; });
}

/*
Slightly changed implementation of Rader's DFT for prime length input. This algorithm computes the adjoint CRT instead of the DFT.
*/
void rader_crt_star_for_primes(complex_vec& output, complex_vec const& input, complex_vec const& precomp_DFT_omega_p, int generator)
{
	int p = input.size() + 1;

	complex_vec permuted_vec = complex_vec(p - 1), zero_padded_vec;

	// the first entry of the output is just the sum of all input entries.
	complex_type out0 = std::accumulate(input.begin(), input.end(), complex_type(REAL_TYPE_ZERO));

	// apply the permutation given by the generator g.
	unitsGeneratorPermutation(permuted_vec, input, generator, p);

	// zero pad to a vector of size equal to the next power of two if necessary.
	if (!isPowerOfTwo(p - 1)){
		zeroPadding(zero_padded_vec, permuted_vec, precomp_DFT_omega_p.size());
	}
	else{
		zero_padded_vec = permuted_vec;
	}

	// compute DFT of this vector via Cooley-Tukey FFT.
	cooley_tukey_fft(zero_padded_vec, zero_padded_vec);

	// component-wise multiply with the pre-computed DFT of the vector of powers of the p-th root of unity.
	std::transform(zero_padded_vec.begin(), zero_padded_vec.end(), precomp_DFT_omega_p.begin(), zero_padded_vec.begin(), std::multiplies<complex_type>());

	// compute IFFT via Cooley-Tukey.
	cooley_tukey_ifft(zero_padded_vec, zero_padded_vec);

	// extract result, i.e., the first p - 1 entries.
	std::copy(zero_padded_vec.begin(), zero_padded_vec.begin() + (p - 1), permuted_vec.begin());

	// undo the permutation.
	unitsGeneratorPermutationInverse(output, permuted_vec, generator, p);

	std::rotate(output.begin(), output.end() - 1, output.end());
	
	output[0] = out0;
}

// Cooley-Tukey FFT for input of length equal to a power of two.
void cooley_tukey_fft(complex_vec& output, complex_vec const& input)
{
	int N = input.size();
	if (!isPowerOfTwo(N)){
		throw "Input length is not a power of 2!";
	}
	else{
		cooley_tukey_fft(output, input, N);
	}
}

// Cooley-Tukey IFFT for input of length equal to a power of two.
void cooley_tukey_ifft(complex_vec& output, complex_vec const& input)
{
	int N = input.size();
	if (!isPowerOfTwo(N)){
		throw "Input length is not a power of 2!";
	}
	else{
		cooley_tukey_ifft(output, input, N);
		// divide result by N
		complex_type N_comp = complex_type(N);
		std::for_each(output.begin(), output.end(), [N_comp](complex_type& e){ e /= N_comp; });
	}
}

// Recursive Cooley-Tukey FFT for input of length equal to a power of two.
void cooley_tukey_fft(complex_vec& output, complex_vec const& input, int N, int start /* = 0*/, int step /* = 1*/)
{
	if (N == 1){ // recursion anchor
		output[start] = input[start];
		return;
	}
	else{
		// recursive step with N / 2.
		cooley_tukey_fft(output, input, N / 2, start, 2 * step); // evenly indexed sub vector
		cooley_tukey_fft(output, input, N / 2, start + step, 2 * step); // oddly indexed sub vector
	}

	// make a temporary working copy of the actually treated sub vector
	complex_vec temp_vec = complex_vec(N);
	int j = 0;
	for (int i = start; i < N*step; i += step)
	{
		temp_vec[j] = output[i];
		j++;
	}

	int n = (N / 2) * step;
	complex_type power = complex_type(REAL_TYPE_ONE), XI = exp(-2 * M_PI * COMP_M_I / real_type(N));
	j = 0;

	// combine both partial results
	for (int i = start; i < n; i += step){
		output[i] = temp_vec[j] + power * temp_vec[j + 1];
		output[i + n] = temp_vec[j] - power * temp_vec[j + 1];
		power *= XI;
		j += 2;
	}
}

// Recursive Cooley-Tukey IFFT for input of length equal to a power of two.
void cooley_tukey_ifft(complex_vec& output, complex_vec const& input, int N, int start /* = 0*/, int step /* = 1*/)
{
	if (N == 1){ // recursion anchor
		output[start] = input[start];
		return;
	}
	else{
		// recursive step with N / 2.
		cooley_tukey_ifft(output, input, N / 2, start, 2 * step); // evenly indexed sub vector
		cooley_tukey_ifft(output, input, N / 2, start + step, 2 * step); // oddly indexed sub vector
	}

	// make a temporary working copy of the actually treated sub vector
	complex_vec temp_vec = complex_vec(N);
	int j = 0;
	for (int i = start; i < N*step; i += step)
	{
		temp_vec[j] = output[i];
		j++;
	}

	int n = (N / 2) * step;
	complex_type power = complex_type(REAL_TYPE_ONE), XI = exp(2 * M_PI * COMP_M_I / real_type(N));
	j = 0;

	// combine both partial results
	for (int i = start; i < n; i += step){
		output[i] = temp_vec[j] + power * temp_vec[j + 1];
		output[i + n] = temp_vec[j] - power * temp_vec[j + 1];
		power *= XI;
		j += 2;
	}
}

// Applies the permutation induced by the generator @param generator of the units of ZZ_p
void unitsGeneratorPermutation(complex_vec& output, complex_vec const& input, int generator, int p)
{
	if (input.size() == p){ // permutes the entries 1,...,p-1 of the input
		int g = 1;
		for (int r = 0; r < p - 1; r++)
		{
			output[r] = input[g];
			g = (g * generator) % p; // traverse through the powers of g mod p.
		}
	}
	else{
		if (input.size() == p - 1){ // permutes the entries 0,...,p-2 of the input
			int g = 1;
			for (int r = 0; r < p - 1; r++)
			{
				output[r] = input[g - 1];
				g = (g * generator) % p; // traverse through the powers of g mod p.
			}
		}
		else{
			throw "Input vector has unexpected size!";
		}
	}
	
}

// Applies the permutation induced by the inverse of the generator @param generator of the units of ZZ_p
void unitsGeneratorPermutationInverse(complex_vec& output, complex_vec const& input, int generator, int p)
{
	if (output.size() == p - 1){ // writes the permuted entries to 0,...,p-2 of the output
		int g_inv = fastModPow(generator, p - 2, p);
		int g = 1;
		for (int r = 0; r < p - 1; r++)
		{
			output[g - 1] = input[r];
			g = (g * g_inv) % p; // traverse through the powers of g^-1 mod p.
		}
	}
	else{
		if (output.size() == p){ // writes the permuted entries to 1,...,p-1 of the output
			int g_inv = fastModPow(generator, p - 2, p);
			int g = 1;
			for (int r = 0; r < p - 1; r++)
			{
				output[g] = input[r];
				g = (g * g_inv) % p; // traverse through the powers of g^-1 mod p.
			}
		}
		else{
			throw "Output vector has unexpected size!";
		}
	}
	
}

void zeroPadding(complex_vec& output, complex_vec const& input, int newSize)
{
	if (newSize < 2 * input.size() - 3){
		throw "The new size is to small.";
	}
	output.resize(newSize);

	// zero padding input
	int M = input.size();
	output[0] = input[0];
	// pad with zeros
	std::fill_n(output.begin() + 1, newSize - M, 0);
	// copy the rest
	std::copy(input.begin() + 1, input.end(), output.begin() + (M + 1));
}

void cyclicPadding(complex_vec& output, complex_vec const& input, int newSize)
{
	if (newSize < 2 * input.size() - 3){
		throw "The new size is to small.";
	}
	output.resize(newSize);

	// cyclic padding input 
	int n = input.size();
	for (int i = 0; i < newSize; i++)
	{
		output[i] = input[i % n];
	}
}



