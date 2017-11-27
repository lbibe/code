

#ifndef PARAMS_HPP
#define PARAMS_HPP

//#include "nfl.hpp"
#include<cstdint>
#include<nfl.hpp>



#define THREADS_NUMBER 2


/**
 * Define the parameters for the polynoms
 */
#define MAX_Q_BITS 130

#define NFL_POLY_COEF_TYPE uint64_t
#define NFL_POLY_N         1024
#define NFL_POLY_Q_BITS    62
using poly_tt = nfl::poly_from_modulus<NFL_POLY_COEF_TYPE, NFL_POLY_N, NFL_POLY_Q_BITS>;
using Poly_t  = nfl::poly_p<typename poly_tt::value_type, poly_tt::degree, poly_tt::nmoduli>;


/**
 * Define the parameters for the Gaussians
 */
using Gauss_t = nfl::gaussian<uint16_t, NFL_POLY_COEF_TYPE, 2>;
using fastGauss_t = nfl::FastGaussianNoise<uint16_t, NFL_POLY_COEF_TYPE, 2>;




#endif