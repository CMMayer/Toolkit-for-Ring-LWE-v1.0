#pragma once

#ifndef DATA_TYPE_SETTINGS_H

#define DATA_TYPE_SETTINGS_H

#include <complex>
#include <vector>

typedef int coordinate_type;	// type for coordinate vectors and integral matrices (int or long)
typedef unsigned int pos_int;

typedef double real_type; // type for representation of the reals. Also for usage in std::complex. 
typedef std::complex<real_type> complex_type;

//vecs
typedef std::vector<coordinate_type> coordinate_vec;
typedef std::vector<real_type> real_vec;
typedef std::vector<complex_type> comp_vec;

#define REAL_TYPE_ZERO real_type(0)
#define REAL_TYPE_ONE real_type(1)

//imaginary unit i
#define COMP_M_I complex_type(REAL_TYPE_ZERO, REAL_TYPE_ONE)

//sqrt(1/2) as a complex number
#define COMP_M_SQRT1_2 complex_type(M_SQRT1_2)

#endif // !DATA_TYPE_SETTINGS_H