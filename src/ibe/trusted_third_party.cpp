#include<cstdint>
#include<cmath>
#include<algorithm>
#include<memory>
#include<thread>
#include<nfl.hpp>

#include "gauss/gaussian.hpp"
#include "params.hpp"
#include "trusted_third_party.hpp"

using namespace Gaussian;

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
struct Ibe::TrustedParty::encapsulation {
	const uint32_t dimension; // n
	const uint32_t param;     // k
	const uint64_t modulus;   // q
	const uint32_t security;  // lambda
	double_t sigma;           // sigma
	uint32_t m; 		  	  // m


	 // Private key: (r, e) --> in NTT
	 Poly_t r[MAX_Q_BITS];
	 Poly_t e[MAX_Q_BITS];

	 // Public key: A, u
	 Poly_t u{nfl::uniform()};
	 Poly_t A[MAX_Q_BITS];  

	// Gaussian
	std::unique_ptr<GaussianSampling> gaussianSampler;
	std::unique_ptr<Gauss_t> gaussianNoise;
	std::unique_ptr<fastGauss_t> fastGaussNoise;

	// SamplePre
	void samplePre(Poly_t * output, const Poly_t * a, const Poly_t & h_1) const noexcept;
};




void Ibe::TrustedParty::encapsulation::samplePre(Poly_t * output, const Poly_t * a, const Poly_t & h_1) const noexcept {
	const Poly_t * p = gaussianSampler->samplePz();

	// Pool of threads
	std::thread pool[THREADS_NUMBER];
	uint32_t thread_shares[THREADS_NUMBER]; 

	// Thread shares
	/*
	const uint32_t thread_share1 = m/THREADS_NUMBER;
	std::fill_n(thread_shares, THREADS_NUMBER, thread_share1);
	for(uint32_t i = 0; i < m - thread_share1*THREADS_NUMBER; ++i) {
		++thread_shares[i];
	}
	*/

	// lambda to compute v
	/*
	auto lambda1 = [](Poly_t * vi, const Poly_t * ai, const Poly_t * pi, const uint32_t length) {
//		vi->ntt_pow_phi();
		for(uint32_t i = 0; i < length; ++i) {
			*vi = (*vi) + ai[i]*pi[i];
		}
	};


	const Poly_t * aIndex = a;
	const Poly_t * pIndex = p;
	Poly_t * vis = new Poly_t[THREADS_NUMBER]{0};
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i] = std::thread(lambda1, vis + i, aIndex, pIndex, thread_shares[i]);
		aIndex += thread_shares[i];
		pIndex += thread_shares[i];
	}
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i].join();
	}

	// compute v
	Poly_t v{0};
//	v.ntt_pow_phi();
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		v = v + *(vis + i);
	}
	v = h_1*(u - v);
	delete[] vis;
	*/


	// compute v
	Poly_t v{0};
	const Poly_t * aIndex = a;
	const Poly_t * pIndex = p;
	for(uint32_t i = 0; i < m; ++i) {
		v = v + (*aIndex)*(*pIndex);
		++aIndex;
		++pIndex;
	}
	v = h_1*(u - v);



	// compute z
	uint64_t vCoefs[dimension];
	uint32_t index = 0;
	for(auto & v_i : v.poly_obj()) {
		vCoefs[index++] = (uint64_t) v_i;
	}

	Poly_t z[param];
	gaussianSampler->sampleGPoly(z, vCoefs);

	// lambda to compute ez, rz and x[2:]
	auto lambda2 = [](Poly_t * outputs, Poly_t * xis, const Poly_t * ei, const Poly_t * ri, const Poly_t * zi, const Poly_t * pi, const uint32_t length) {
//		outputs[0].ntt_pow_phi();
//		outputs[1].ntt_pow_phi();
		for(uint32_t i = 0; i < length; ++i) {
			xis[i]     = pi[i] + zi[i];
			outputs[0] = outputs[0] + ei[i]*zi[i];
			outputs[1] = outputs[1] + ri[i]*zi[i];
		}
	};

	// Thread shares
	const uint32_t thread_share2 = param/THREADS_NUMBER;
	std::fill_n(thread_shares, THREADS_NUMBER, thread_share2);
	for(uint32_t i = 0; i < param - thread_share2*THREADS_NUMBER; ++i) {
		++thread_shares[i];
	}


	pIndex          = p + 2;
	Poly_t * xIndex = output + 2;
	const Poly_t * eIndex = e;
	const Poly_t * rIndex = r;
	const Poly_t * zIndex = z;
	Poly_t * erz = new Poly_t[THREADS_NUMBER*2]{0};
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		const uint32_t share = thread_shares[i];
		pool[i] = std::thread(lambda2, erz + i*2, xIndex, eIndex, rIndex, zIndex, pIndex, share);
		xIndex += share;
		eIndex += share;
		rIndex += share;
		zIndex += share;
		pIndex += share;
	}
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i].join();
	}


	// compute ez, rz and x[2:]
	Poly_t ez{0};
	Poly_t rz{0};
//	ez.ntt_pow_phi();
//	rz.ntt_pow_phi();
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		ez = ez + erz[i*2];
		rz = rz + erz[1 + i*2];
	}
	delete[] erz;

	// compute x
	output[0] = p[0] + ez;
	output[1] = p[1] + rz;

	delete[] p;
}







/* --------------------------------------------------------------------------------------------------- */
/* Public Interface ---------------------------------------------------------------------------------- */
/* @Override */
Ibe::TrustedParty::TrustedParty(const uint32_t n, const uint32_t k, const uint64_t q, const double_t sigma, const uint32_t lambda) noexcept : \
impl(new encapsulation {.dimension = n, .param = k, .modulus = q, .security = lambda}) {
	impl->sigma = sigma;
	impl->m = k + 2;


	// prebuild the Gaussian noise
	impl->fastGaussNoise.reset(new fastGauss_t(sigma, lambda, n));
	impl->gaussianNoise.reset(new Gauss_t(impl->fastGaussNoise.get()));

	// For the Public key
	impl->u.ntt_pow_phi();
	impl->A[0].set(1);
	impl->A[0].ntt_pow_phi();
}

/* @Override */
Ibe::TrustedParty::~TrustedParty(void) noexcept {

}


/* @Override */
void Ibe::TrustedParty::generateMasterKey(void) noexcept {
	// Get k
	const uint32_t k = impl->param;

	// compute a
	Poly_t a{nfl::uniform()};
	a.ntt_pow_phi();

	// Pool of threads
	std::thread pool[THREADS_NUMBER];

	// Thread shares
	uint32_t thread_shares[THREADS_NUMBER]; 
	const uint32_t thread_share = (k*2)/THREADS_NUMBER;
	std::fill_n(thread_shares, THREADS_NUMBER, thread_share);
	for(uint32_t i = 0; i < k*2 - thread_share*THREADS_NUMBER; ++i) {
		++thread_shares[i];
	}



	// Lambda for ^r and ê
	auto lambda1 = [this](Poly_t * output, const uint32_t length) {
		Gauss_t noise = *(impl->gaussianNoise.get());
		for(uint32_t i = 0; i < length; ++i) {
			output->set(noise);
			output->ntt_pow_phi();
			++output;
		}
	};

	// compute ^r
	Poly_t * rIndex = impl->r;
	for(uint32_t i = 0; i < THREADS_NUMBER/2; ++i) {
		pool[i] = std::thread(lambda1, rIndex, thread_shares[i]);
		rIndex += thread_shares[i];
	}
	/*
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i].join();
	}
	*/

	// compute ê
	Poly_t * eIndex = impl->e;
	for(uint32_t i = THREADS_NUMBER/2; i < THREADS_NUMBER; ++i) {
		pool[i] = std::thread(lambda1, eIndex, thread_shares[i]);
		eIndex += thread_shares[i];
	}
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i].join();
	}

	// Lambda for A
	auto lambda2 = [&a](Poly_t * output, const Poly_t * r, const Poly_t * e, const uint32_t length) {
		for(uint32_t i = 0; i < length; ++i) {
			*output = - (a*(*r) + (*e));
			++output;
			++r;
			++e;		
		}
	};


	// compute A
	rIndex = impl->r;
	eIndex = impl->e;
	Poly_t * AIndex = (impl->A + 2);
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		const uint32_t share = thread_shares[i];
		pool[i]  = std::thread(lambda2, AIndex, rIndex, eIndex, share);
		AIndex  += share;
		rIndex  += share;
		eIndex  += share;
	}
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i].join();
	}

	impl->A[1] = std::move(a);
}



/* @Override */
Poly_t * Ibe::TrustedParty::getPublicKey(Poly_t * pk) noexcept {
	memcpy(pk, &(impl->u), sizeof(Poly_t));
	return impl->A;
}

/* @Override */
void Ibe::TrustedParty::setGaussian(const double_t s) noexcept {
	impl->gaussianSampler.reset(new GaussianSampling(impl->dimension, impl->param, impl->modulus, impl->sigma, s, 2*impl->sigma, impl->security, impl->r, impl->e));
	std::cout << "End of gaussian" << std::endl;
}


/* @Override */
void Ibe::TrustedParty::preCompute(const uint32_t n) noexcept {
	impl->gaussianSampler->preCompute(n);
}

/* @Override */
void Ibe::TrustedParty::extract(Poly_t * output, const Poly_t * a, const Poly_t & h_1) const noexcept {
	impl->samplePre(output, a, h_1);
}

/* Public Interface End ------------------------------------------------------------------------------ */
/* --------------------------------------------------------------------------------------------------- */