
#ifndef GAUSSIAN_HPP
#define GAUSSIAN_HPP

#include<cstdint>
#include<cmath>
#include<memory>
#include<algorithm>
#include<thread>
#include "params.hpp"


namespace Gaussian {
	class GaussianSampling {

		public:

			/** 
			 * @Default Constructor
			 */
			GaussianSampling(void) noexcept;

			/** 
			 * @Constructor
			 */
			GaussianSampling(const uint32_t n, const uint32_t k, const uint64_t q, const double_t s, const double_t ps, const double_t alpha, const uint32_t lambda, const Poly_t * r, const Poly_t * e) noexcept;

			/** 
			 * @Move Constructor
			 */
			GaussianSampling(GaussianSampling && other) noexcept;

			/** 
			 * @Destructor
			 */
			~GaussianSampling(void) noexcept;


			/** 
			 * @Move assignment
			 */
			GaussianSampling & operator=(GaussianSampling && other) noexcept;


			/**
			 *
			 */
			void preCompute(const uint32_t n) noexcept;


			/**
			 *
			 */
			void sampleGPoly(Poly_t * polys, const uint64_t * coset) const noexcept;

			/**
			 *
			 */
			Poly_t * samplePz(void) const noexcept;


			/**
			 * @test
			 */
//			void test(nfl::poly_p<NFL_POLY_COEF_TYPE, NFL_POLY_N, NFL_POLY_Q_BITS> * polys) noexcept;

		private:
			struct encapsulation;
			std::unique_ptr<encapsulation> impl;

			void G_sampleGsMonoThread(int32_t * output, const uint64_t * coset, const uint32_t length) const noexcept;
	};
}

#endif
