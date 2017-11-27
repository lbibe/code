#ifndef USER_HPP
#define USER_HPP

#include<cstdint>
#include<cmath>
#include<memory>
#include "params.hpp"
#include "trusted_third_party.hpp"


namespace Ibe {

	class User {
		public:

			/** 
			 * @Constructor
			 */
			User(const uint32_t n, const uint32_t k, const uint64_t q, const double_t tau, const double_t gamma, const uint32_t lambda, \
				 const NFL_POLY_COEF_TYPE id, TrustedParty * trustedParty) noexcept;

			/** 
			 * @Destructor
			 */
			~User(void) noexcept;


			/**
			 * @
			 */
			void extractPrivateKey(void) const noexcept;


			/**
			 * @
			 */
			void setEncrypt(const NFL_POLY_COEF_TYPE targetId) const noexcept;


			void encrypt(Poly_t * output, const Poly_t & msg) const noexcept;


			void decrypt(Poly_t * output, const Poly_t * cipher) const noexcept;


		private:
			struct encapsulation;
			std::unique_ptr<encapsulation> impl;
	};

}


#endif