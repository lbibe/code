#include<iostream>
#include <chrono>
#include<thread>
#include<algorithm>
#include<complex>
#include "util/fast_fft.hpp"
#include "ibe/user.hpp"
#include "ibe/trusted_third_party.hpp"

#include<nfl.hpp>

#define sigma_1024_ibe 5
#define zeta_1024_ibe 6360
#define tau_1024_ibe 5


using namespace Ibe;
int main(void) {
	// Pool of threads
	std::thread pool[2];


	mpz_t mod;
    mpz_init(mod);
    mpz_set(mod, Poly_t::moduli_product());
    const uint32_t modulo = mpz_get_ui(mod);

	std::cout << "Before TrustedParty" << std::endl;
	TrustedParty thirdParty{NFL_POLY_N, NFL_POLY_Q_BITS, modulo, sigma_1024_ibe, 80};
	thirdParty.generateMasterKey();
	thirdParty.setGaussian(zeta_1024_ibe);
	thirdParty.preCompute(2);


	std::cout << "Before User" << std::endl;
	User sender{NFL_POLY_N, NFL_POLY_Q_BITS, modulo, tau_1024_ibe, 147.47, 80, 17, &thirdParty};
	sender.extractPrivateKey();
	sender.setEncrypt(505);


	Poly_t m = nfl::uniform();
	for (auto &m_i : m.poly_obj()) {
		m_i %= 2;
	}
	Poly_t * ciphertext = new Poly_t[(NFL_POLY_Q_BITS + 3)*2];
	sender.encrypt(ciphertext, m);


	auto lambdaEnc = [&sender, &m](Poly_t * ciphertext, const uint32_t length) {
		for(uint32_t i = 0; i < length; ++i) {
			sender.encrypt(ciphertext, m);
		}
	};
	pool[0] = std::thread(lambdaEnc, ciphertext, 100);
	pool[1] = std::thread(lambdaEnc, ciphertext + NFL_POLY_Q_BITS + 3, 100);
	auto start = std::chrono::high_resolution_clock::now();
	pool[0].join();
	pool[1].join();
	auto end = std::chrono::high_resolution_clock::now();
	auto timing = std::chrono::duration<double, std::milli>(end - start);
	std::cout << "Encryption Time of 200 messages: " << timing.count() << " ms " << std::endl;


	Poly_t * decrypted = new Poly_t[2];
	sender.decrypt(decrypted, ciphertext);

	auto lambdaDec = [&sender, ciphertext](Poly_t * output, const uint32_t length) {
		for(uint32_t i = 0; i < length; ++i) {
			sender.decrypt(output, ciphertext);
		}
	};

	pool[0] = std::thread(lambdaDec, decrypted, 1000);
	pool[1] = std::thread(lambdaDec, decrypted + 1, 1000);
	start = std::chrono::high_resolution_clock::now();
	pool[0].join();
	pool[1].join();
	end = std::chrono::high_resolution_clock::now();
	timing = std::chrono::duration<double, std::milli>(end - start);
	std::cout << "Decryption Time of 2000 ciphertexts: " << timing.count() << " ms " << std::endl;

	return 0;
}