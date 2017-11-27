
#ifndef TRUSTED_THIRD_PARTY_HPP
#define TRUSTED_THIRD_PARTY_HPP

#include<cstdint>
#include<cmath>
#include<memory>
#include "params.hpp"



namespace Ibe {

	class TrustedParty {

		public:

			/** 
			 * @Constructor
			 */
			TrustedParty(const uint32_t n, const uint32_t k, const uint64_t q, const double_t sigma, const uint32_t labmda) noexcept;

			/** 
			 * @Destructor
			 */
			~TrustedParty(void) noexcept;



			/** 
			 * @
			 */
			void generateMasterKey(void) noexcept;


			/** 
			 * @
			 */
			Poly_t * getPublicKey(Poly_t * pk) noexcept;


			/** 
			 * @
			 */
			void setGaussian(const double_t s) noexcept;


			/** 
			 * @
			 */
			void preCompute(const uint32_t n) noexcept;


			/** 
			 * @
			 */
			void extract(Poly_t * output, const Poly_t * a, const Poly_t & h_1) const noexcept;


		private:
			struct encapsulation;
			std::unique_ptr<encapsulation> impl;
	};
}





#endif