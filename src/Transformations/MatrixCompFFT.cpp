#include <algorithm>

#include "Transformations/MatrixCompFFT.h"
#include "Math util.h"
#include "Efficient DFT.h"

typedef MatrixCompFFT::MatrixType trans_type;

MatrixCompFFT::MatrixCompFFT(trans_type matrixType, int p) :
AbstractVectorTransformation(p),
trans_(matrixType)
{
	generator = findGeneratorOfZZpUnits(p);
	int generator_inv = fastModPow(generator, p - 2, p);
	complex_type rootOfUnity = computeRootOfUnity(p);

	complex_vec root_powers = complex_vec(p - 1);
	switch (trans_)
	{
	case trans_type::CRT_P:
		dim_ = p - 1;
	case trans_type::DFT_P:
		for (int r = 0; r < p - 1; r++)
		{
			// traverse through the roots of unity according to the order 
			// induced by the inverse of the generator
			int exp = fastModPow(generator_inv, r, p);
			root_powers[r] = pow(rootOfUnity, exp);
		}
		break;

	case trans_type::DFT_P_INV:
		for (int r = 0; r < p - 1; r++)
		{
			// traverse through the roots of unity according to the order 
			// induced by the inverse of the generator
			int exp = fastModPow(generator_inv, r, p);
			root_powers[r] = pow(rootOfUnity, -exp); // take inverse of the root
		}
		break;

	case trans_type::CRT_P_STAR:
		dim_ = p - 1;
	case trans_type::DFT_P_STAR:
		for (int r = 0; r < p - 1; r++)
		{
			// traverse through the roots of unity according to the order 
			// induced by the inverse of the generator
			int exp = fastModPow(generator_inv, r, p);
			root_powers[r] = conj(pow(rootOfUnity, exp)); // take conjugate of the root
		}
		break;


	case trans_type::DFT_P_STAR_INV:
		for (int r = 0; r < p - 1; r++)
		{
			// traverse through the roots of unity according to the order 
			// induced by the inverse of the generator
			int exp = fastModPow(generator_inv, r, p);
			// take conjugate of the inverse, which is the original root, since conjugation 
			// and inversion coincide for roots of unity
			root_powers[r] = pow(rootOfUnity, exp); 
		}
		break;

	default:
		break;
	}

	precomp_DFT_omega_p_ = std::make_unique<complex_vec>();
	if (!isPowerOfTwo(p - 1)){ // cyclic padding if size is not a power of two
		int power_of_two = getNextPowerOfTwo(2 * (p - 1) - 3);
		precomp_DFT_omega_p_->resize(power_of_two);
		cyclicPadding(*precomp_DFT_omega_p_, root_powers, power_of_two);
	}
	else{
		*precomp_DFT_omega_p_ = root_powers;
	}
	// pre-compute FFT of the roots.
	cooley_tukey_fft(*precomp_DFT_omega_p_, *precomp_DFT_omega_p_);
}

MatrixCompFFT::~MatrixCompFFT(){}

// Apply this matrix to @param vec and store the result in x.
void MatrixCompFFT::applyToVector(entry_type_vec& x, entry_type_vec const& vec) const
{
	switch (trans_)
	{
	case trans_type::CRT_P:
		rader_crt_for_primes(x, vec, *precomp_DFT_omega_p_, generator);
		break;

	case trans_type::CRT_P_STAR:
		rader_crt_star_for_primes(x, vec, *precomp_DFT_omega_p_, generator);
		break;

	case trans_type::DFT_P:
	case trans_type::DFT_P_STAR: // the vector precomp_DFT_omega_p_ makes sure
		// that Rader DFT corresponds to an application of DFT_p or DFT_p^*
		rader_dft_for_primes(x, vec, *precomp_DFT_omega_p_, generator);
		break;

	case trans_type::DFT_P_INV:
	case trans_type::DFT_P_STAR_INV:
	{
		// the vector precomp_DFT_omega_p_ makes sure
		// that Rader DFT corresponds to an application of p*DFT_p^-1 or p*(DFT_p^*)^-1
		rader_dft_for_primes(x, vec, *precomp_DFT_omega_p_, generator);
		entry_type scale = entry_type(dim_);
		std::for_each(x.begin(), x.end(), [scale](entry_type& entry){ entry /= scale; });
		break;
	}

	default:
		break;
	}
}
