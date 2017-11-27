
#include<cstdint>
#include<cmath>
#include<random>
#include<algorithm>
#include<memory>
#include<thread>
#include<nfl.hpp>

#include "user.hpp"
#include "params.hpp"



/* --------------------------------------------------------------------------------------------------- */
/* Const Values -------------------------------------------------------------------------------------- */

/* Const Values End ---------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------------------------- */



/* --------------------------------------------------------------------------------------------------- */
/* Static Functions ---------------------------------------------------------------------------------- */

/* Static Functions End ------------------------------------------------------------------------------ */
/* --------------------------------------------------------------------------------------------------- */



/*
 * Private fields
 */
 struct Ibe::User::encapsulation {
 	//
 	const uint32_t dimension; // n
	const uint32_t param;     // k
	const uint64_t modulus;   // q
	uint32_t m; 		  	  // m

 	// ID
 	NFL_POLY_COEF_TYPE id;
 	Poly_t h_1{0};
 	Poly_t a_id[MAX_Q_BITS];

 	// encryption status
 	mutable Poly_t encryption_a_id[MAX_Q_BITS];


 	// trusted party & public/private key
 	Poly_t * mpkA;
 	Poly_t mpkU{0};
 	Poly_t sk[MAX_Q_BITS];
 	TrustedParty * trustedParty;

 	// Related to Trapdoor
 	Poly_t gs[MAX_Q_BITS];

 	// Gaussian errors
 	std::unique_ptr<Gauss_t> gaussianNoiseE0;
 	std::unique_ptr<Gauss_t> gaussianNoiseE1;
 	std::unique_ptr<Gauss_t> gaussianNoiseE2;
	std::unique_ptr<fastGauss_t> fastGaussNoiseE0;
	std::unique_ptr<fastGauss_t> fastGaussNoiseE1;
	std::unique_ptr<fastGauss_t> fastGaussNoiseE2;

	// 
	std::independent_bits_engine<std::mt19937_64, NFL_POLY_Q_BITS, NFL_POLY_COEF_TYPE> randomGen;

	// Functions
	void hash(Poly_t * output, const NFL_POLY_COEF_TYPE x);
	void compute_a_id(Poly_t * output, const Poly_t & h);
	NFL_POLY_COEF_TYPE inverse(const NFL_POLY_COEF_TYPE a);
	void inverse_poly(Poly_t * output, Poly_t & input);
 };



void Ibe::User::encapsulation::hash(Poly_t * output, const NFL_POLY_COEF_TYPE x) {	
	randomGen.seed(x);
	for(auto & output_i : output->poly_obj()) {
		output_i = randomGen();
	}
}


void Ibe::User::encapsulation::compute_a_id(Poly_t * output, const Poly_t & h) {	
//	output[0] = mpkA[0]; output[1] = mpkA[1];
	memcpy(output, mpkA, sizeof(Poly_t)*2);

	// Pool of threads
	std::thread pool[THREADS_NUMBER];

	// Thread shares
	uint32_t thread_shares[THREADS_NUMBER]; 
	const uint32_t thread_share = param/THREADS_NUMBER;
	std::fill_n(thread_shares, THREADS_NUMBER, thread_share);
	for(uint32_t i = 0; i < param - thread_share*THREADS_NUMBER; ++i) {
		++thread_shares[i];
	}

	// Lambda for a_id
	auto lambda = [&h](Poly_t * outputs, const Poly_t * gs, const Poly_t * a, const uint32_t length) {
		for(uint32_t i = 0; i < length; ++i) {
			*outputs = h*(*gs) + (*a);
			++outputs;
			++gs;
			++a;		
		}
	};

	// compute a_id
	const Poly_t * gsIndex = gs;
	const Poly_t * AIndex = (mpkA + 2);
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		const uint32_t share = thread_shares[i];
		pool[i]  = std::thread(lambda, output, gsIndex, AIndex, share);
		output  += share;
		AIndex  += share;
		gsIndex += share;
	}
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i].join();
	}
}


NFL_POLY_COEF_TYPE Ibe::User::encapsulation::inverse(const NFL_POLY_COEF_TYPE a) {
    NFL_POLY_COEF_TYPE u0 = 1, v0 = 0;
    NFL_POLY_COEF_TYPE u1 = 0, v1 = 1;
    NFL_POLY_COEF_TYPE res0 = a, res1 = modulus, quo = 0;
    NFL_POLY_COEF_TYPE tmp;

    while (res1 != 0) {
      quo = res0 / res1;
      tmp = res0, res0 = res1, res1 = tmp % res1;
      tmp = u0, u0 = u1, u1 = tmp - quo*u1;
      tmp = v0, v0 = v1, v1 = tmp - quo*u1;
    }
    return (u0 > 0) ? u0 : u0 + modulus;
}

void Ibe::User::encapsulation::inverse_poly(Poly_t * output, Poly_t & input) {
	/*
	auto inputs = input.poly2mpz();
	const uint32_t n = inputs.size();

	NFL_POLY_COEF_TYPE invs[n];
	for(uint32_t i = 0; i < n; ++i) {
		invs[i] = inverse(mpz_get_ui(inputs[i]));
	}
    output->set(invs, invs + n);
    */
    uint32_t i = 0;
    NFL_POLY_COEF_TYPE invs[dimension];
	for(auto & input_i : input.poly_obj()) {
		invs[i++] = inverse(input_i);
	}
	output->set(invs, invs + dimension, false);
}


 /* --------------------------------------------------------------------------------------------------- */
/* Public Interface ---------------------------------------------------------------------------------- */
/* @Override */
Ibe::User::User(const uint32_t n, const uint32_t k, const uint64_t q, const double_t tau, const double_t gamma, const uint32_t lambda, \
	            const NFL_POLY_COEF_TYPE id, TrustedParty * trustedParty) noexcept : \
impl(new encapsulation {.dimension = n, .param = k, .modulus = q}) {
	// initialize global parameters
	impl->id = id;
	impl->m = k + 2;

	// initialize trusted party & public keys
	impl->mpkA = trustedParty->getPublicKey(&impl->mpkU);
	impl->trustedParty = trustedParty;


	// prebuild the Gaussian noises
	impl->fastGaussNoiseE0.reset(new fastGauss_t(tau, lambda, n));
	impl->gaussianNoiseE0.reset(new Gauss_t(impl->fastGaussNoiseE0.get()));

	impl->fastGaussNoiseE1.reset(new fastGauss_t(gamma, lambda, n));
	impl->gaussianNoiseE1.reset(new Gauss_t(impl->fastGaussNoiseE1.get()));

	impl->fastGaussNoiseE2.reset(new fastGauss_t(tau, lambda, n));
	impl->gaussianNoiseE2.reset(new Gauss_t(impl->fastGaussNoiseE2.get()));


	// build gs and NTT
	NFL_POLY_COEF_TYPE gi = 1;
	for(uint32_t i = 0; i < k; ++i) {
		impl->gs[i].set(gi);
		impl->gs[i].ntt_pow_phi();
		gi <<= 1;
	}


	// compute h
	Poly_t h{0};
	impl->hash(&h, id);


	// inverse h
	impl->inverse_poly(&(impl->h_1), h);


	// compute a_id
	impl->compute_a_id(impl->a_id, h);

}

/* @Override */
Ibe::User::~User(void) noexcept {

}


/* @override */
void Ibe::User::extractPrivateKey(void) const noexcept {
	impl->trustedParty->extract(impl->sk, impl->a_id, impl->h_1);
}


/* @override */
void Ibe::User::setEncrypt(const NFL_POLY_COEF_TYPE targetId) const noexcept {
	// compute h
	Poly_t h{0};
	impl->hash(&h, targetId);

	// compute a_id
	impl->compute_a_id(impl->encryption_a_id, h);
}

/* @override */
void Ibe::User::encrypt(Poly_t * output, const Poly_t & msg) const noexcept {
	// prepare the plaintext
	Poly_t plaintext{msg};
	for (auto & plaintext_i : plaintext.poly_obj()) {
		plaintext_i *= (impl->modulus/2);
    }
    plaintext.ntt_pow_phi();

    // uniform noise
    Poly_t s = nfl::uniform();
    s.ntt_pow_phi();
	
    // gaussian noise e0
    const Poly_t * A_idIndex = impl->encryption_a_id;
    Gauss_t noiseE0 = *(impl->gaussianNoiseE0.get());
    for(uint32_t i = 0; i < impl->m - impl->param; ++i) {
    	output->set(noiseE0);
    	output->ntt_pow_phi();
    	*output = *output + s*(*A_idIndex);
    	++output;
    	++A_idIndex;
    }

    // gaussian noise e1
    Gauss_t noiseE1 = *(impl->gaussianNoiseE1.get());
    for(uint32_t i = 0; i < impl->param; ++i) {
    	output->set(noiseE1);
    	output->ntt_pow_phi();
    	*output = *output + s*(*A_idIndex);
    	++output;
    	++A_idIndex;
    }

    // Pool of threads
    /*
    std::thread pool[THREADS_NUMBER];

    // Thread shares
	uint32_t thread_shares[THREADS_NUMBER]; 
	const uint32_t thread_share = impl->param/THREADS_NUMBER;
	std::fill_n(thread_shares, THREADS_NUMBER, thread_share);
	for(uint32_t i = 0; i < impl->param - thread_share*THREADS_NUMBER; ++i) {
		++thread_shares[i];
	}

    // Lambda for gaussian noise e1
	auto lambda = [this,&s](Poly_t * output, const Poly_t * a_id, const uint32_t length) {
		Gauss_t noiseE1 = *(impl->gaussianNoiseE1.get());
		for(uint32_t i = 0; i < length; ++i) {
			output->set(noiseE1);
			output->ntt_pow_phi();
			*output = *output + s*(*a_id);
			++output;
			++a_id;
		}
	};


	// gaussian noise e1
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i] = std::thread(lambda, output, A_idIndex, thread_shares[i]);
		output    += thread_shares[i];
		A_idIndex += thread_shares[i];
	}
	*/

    // final
    Gauss_t noiseE2 = *(impl->gaussianNoiseE2.get());
	output->set(noiseE2);
	output->ntt_pow_phi();
	*output = *output + impl->mpkU*s + plaintext;

	/*
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i].join();
	}
	*/
}


/* @override */
void Ibe::User::decrypt(Poly_t * output, const Poly_t * cipher) const noexcept {
	Poly_t bx{0};

	const Poly_t * sk_id = impl->sk;
	for(uint32_t i = 0; i < impl->m; ++i) {
		bx = bx + (*sk_id)*(*cipher);
		++sk_id;
		++cipher;
	} 


	*output = *cipher - bx;
	output->invntt_pow_invphi();

	/*
	const uint32_t modulus_4  = impl->modulus >> 2;
	const uint32_t modulus_34 = 3*modulus_4;
	for(auto & output_i : output->poly_obj()) {
		output_i = ((output_i <= modulus_4) || (output_i >= modulus_34)) ? 0 : 1;
	}
	*/
	const uint32_t modulus_4  = impl->modulus >> 2;
	const uint32_t modulus_34 = 3*modulus_4;
	for(auto & output_i : output->poly_obj()) {
		output_i = ((output_i >= modulus_4) && output_i <= modulus_34) ? 1 : 0;
	}
}


 /* Public Interface End ------------------------------------------------------------------------------ */
/* --------------------------------------------------------------------------------------------------- */