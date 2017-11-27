
#include<cstdint>
#include<cstring>
#include<memory>
#include<cmath>
#include<random>
#include<stack>
#include<algorithm>
#include<thread>
#include<mutex>
#include<atomic>
#include "gaussian.hpp"
#include "params.hpp"
#include "util/fast_fft.hpp"


using namespace FFT;


/* --------------------------------------------------------------------------------------------------- */
/* Const Values -------------------------------------------------------------------------------------- */

/* Const Values End ---------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------------------------------- */
/* Static Functions ---------------------------------------------------------------------------------- */
/*
static void merge(double_t * output, const double_t * input_1, const double_t * input_2, const uint32_t length) {
	for(uint32_t i = 0, j = 0; i < length; ++i) {
		output[j++] = input_2[i];
		output[j++] = input_1[i];
	}
}

static void split(double_t * output_1, double_t * output_2, const double_t * input, const uint32_t length) {
	for(uint32_t i = 0, j = 0; i < length; ++j) {
		output_2[j] = input[i++]; 
		output_1[j] = input[i++];
	}
}
*/

/*
TO DELETE
static Poly_t conjugate(Poly_t & polyIn) {
	auto coefsIn = polyIn.poly2mpz();

	std::array<mpz_t, NFL_POLY_N> coefsOut{};
	mpz_set(coefsOut[0], coefsIn[0]);

	const uint32_t n = coefsIn.size();
	for(uint32_t i = 1; i < n; ++i) {
		mpz_neg(coefsOut[i], coefsIn[n - i]);
	}

	Poly_t polyOut;
	polyOut.mpz2poly(coefsOut);
	return polyOut;
}
*/


inline static void getBits(uint8_t * bits, uint64_t n, uint32_t length) {
	uint32_t mask = 1;
	for(uint32_t i = 0; i < length; ++i) {
		bits[i] = ((n & mask)== 0) ? 0 : 1;
		mask <<= 1;
	}
}

/* Static Functions End ------------------------------------------------------------------------------ */
/* --------------------------------------------------------------------------------------------------- */


/*
 * Private fields
 */
struct Gaussian::GaussianSampling::encapsulation {
	// Global variables
	const uint32_t dimension;      // n
	const uint32_t param;          // k
	uint8_t modulus[MAX_Q_BITS];   // q
	double_t sigma;                // s

	// utility for fft
	const FastFFT<NFL_POLY_N> fastFft{};

	std::mt19937_64 * randomGenerator;
	std::unique_ptr<std::mt19937_64> randomGenerator_t;
	int64_t sampleZ(const double_t sigma, const double_t center) const noexcept;

	// Variables associated to sampleG
//	const int32_t G_b = 2;
	double_t G_l[MAX_Q_BITS];
	double_t G_h[MAX_Q_BITS];
	double_t G_d[MAX_Q_BITS];
	mutable std::stack<int32_t*> G_precomputedPerturbations;
	mutable std::mutex G_precomputed_mutex;

	// Functions associated to sampleG
	void G_preCompute(void) const noexcept;
	int32_t * G_Perturb(void) const noexcept;
	inline void G_sampleD(int32_t * output, const double_t sigma, const double_t * c) const noexcept;
	void G_sampleG(int32_t * output, const uint64_t coset_i) const noexcept;


	// Variables associated to sampleP
	double_t P_z0;                // z
	double_t P_z1;                // sqrt(s^2 - alpha^2)
	double_t P_z2;                // (-alpha^2)/(s^2 - alpha^2)
	const Poly_t * P_r;                 // ^r
	const Poly_t * P_e;                 // ê
	std::complex<double_t> P_a[NFL_POLY_N];     // a
	std::complex<double_t> P_b[NFL_POLY_N];     // b
	std::complex<double_t> P_d[NFL_POLY_N];     // d
	std::unique_ptr<Gauss_t> P_gaussianNoise;
	std::unique_ptr<fastGauss_t> P_fastGaussNoise;

	mutable std::stack<Poly_t *> P_precomputedPerturbations;
	mutable std::mutex P_precomputed_mutex;


	// Functions associated to sampleP
	void P_preCompute(void) const noexcept;
	void P_samplePz(Poly_t * output) const noexcept;
	void P_sampleFz(std::complex<double_t> * output, const std::complex<double_t> * f, const std::complex<double_t> * c, const uint32_t length) const noexcept;
	void P_sample2z(std::complex<double_t> * output_1, std::complex<double_t> * output_2, const std::complex<double_t> * a, \
					 const std::complex<double_t> * b, const std::complex<double_t> * d, const std::complex<double_t> * c, const uint32_t length) const noexcept;
};


int64_t Gaussian::GaussianSampling::encapsulation::sampleZ(const double_t sigma, const double_t center) const noexcept {
	// normal distribution
	std::normal_distribution<> normalDist(center, sigma);

	// sample
	return (int64_t) round(normalDist(*randomGenerator));
}


int32_t * Gaussian::GaussianSampling::encapsulation::G_Perturb(void) const noexcept {
	std::lock_guard<std::mutex> lock(G_precomputed_mutex);
	int32_t * result = G_precomputedPerturbations.top();
	G_precomputedPerturbations.pop();
	return result;
}

void Gaussian::GaussianSampling::encapsulation::G_preCompute(void) const noexcept {
	int32_t * perturbation{new int32_t[param]};
	
	// perturb(sigma) in SampleG
	// -----------------------------------------------
	double_t sigmas[param];
	for(uint32_t i = 0; i < param; ++i) {
		sigmas[i] = sigma/G_l[i];
	}

	double_t beta = 0.0;
	int32_t z[param];
	for(uint32_t i = 0; i < param; ++i) {
		double_t c  = beta/G_l[i];
		int32_t  cz = (int32_t) floor(c);
		z[i] = cz + sampleZ(sigmas[i], c - cz);
		beta = -z[i]*G_h[i];
	}

	perturbation[0] = z[0] + (z[0]<<2) + (z[1]<<1);
	for(uint32_t i = 1; i < param; ++i) {
		perturbation[i] = (z[i-1] + (z[i]<<1) + z[i+1])<<1;
	}
	// -----------------------------------------------
	// add to the precomputed perturbations
	std::lock_guard<std::mutex> lock(G_precomputed_mutex);
	G_precomputedPerturbations.push(perturbation);
}

inline void Gaussian::GaussianSampling::encapsulation::G_sampleD(int32_t * output, const double_t sigmaa, const double_t * c) const noexcept {
	const double_t c_d = -c[param - 1]/G_d[param - 1];
	const int32_t c_dz = (int32_t) floor(c_d);
	output[param - 1] = c_dz + sampleZ(sigmaa/G_d[param - 1], c_d - c_dz);


	double_t cs[param];
	const int32_t zk_1 = output[param - 1];
	for(uint32_t i = 0; i < param; ++i) {
		cs[i] = zk_1*G_d[i] - c[i]; 
	}

	for(uint32_t i = 0; i < param - 1; ++i) {
		const int32_t cs_iz = (int32_t) floor(cs[i]); 
		output[i] = cs_iz + sampleZ(sigmaa, cs[i] - cs_iz); 
	}
}


void Gaussian::GaussianSampling::encapsulation::G_sampleG(int32_t * output, const uint64_t coset_i) const noexcept {
	// get bits of coset_i
	uint8_t u[param];
	getBits(u, coset_i, param);

	// get perturbations
	int32_t * p = G_Perturb();

	// compute cs
	double_t cs[param];
	cs[0] = (u[0] - p[0])/2;
	for(uint32_t i = 1; i < param; ++i) {
		cs[i] = (cs[i - 1] + u[i] - p[i])/2;
	}

	// compute z
	int32_t z[param];
	G_sampleD(z, sigma, cs);

	// compute t
	output[0] = (z[0]<<1) + u[0];
	for(uint32_t i = 1; i < param -1; ++i) {
		output[i] = (z[i]<<1) - z[i-1] + modulus[i]*z[i-1] + u[i];
	}
	output[param-1] = modulus[param-1]*z[param-1] - z[param-2] + u[param-1];

	delete[] p;
}



void Gaussian::GaussianSampling::encapsulation::P_sampleFz(std::complex<double_t> * output, const std::complex<double_t> * f, \
														 const std::complex<double_t> * c, const uint32_t length) const noexcept {
	// if dim(f) = 1 return SampleZ(f, c)
	if (length == 1) {
		*output = (double_t) sampleZ(f[0].real(), c[0].real());
	}
	else {
		const uint32_t length_2 = length >> 1;
		// let f(x) = f0(x^2) + x.f1(x^2)
		std::complex<double_t> f0[length_2];
		std::complex<double_t> f1[length_2];
		fastFft.ffsplit(f0, f1, f, length);

		// (q0, q1) Sample2Z(f0, f1, f0, c)
		std::complex<double_t> * q0 = new std::complex<double_t>[length_2];
		std::complex<double_t> * q1 = new std::complex<double_t>[length_2];
		P_sample2z(q0, q1, f0, f1, f0, c, length_2);

		// let q(x) = q0(x^2) + x.q1(x^2)
		fastFft.ffmerge(output, q0, q1, length_2);
		delete[] q0;
		delete[] q1;
	}


}


void Gaussian::GaussianSampling::encapsulation::P_sample2z(std::complex<double_t> * output_1, std::complex<double_t> * output_2, const std::complex<double_t> * a, \
			  const std::complex<double_t> * b, const std::complex<double_t> * d, const std::complex<double_t> * c, const uint32_t length) const noexcept {

	// let c(x) = c0(x^2) + x.c1(x^2)
	std::complex<double_t> c0[length];
	std::complex<double_t> c1[length];
	fastFft.ffsplit(c0, c1, c, length << 1);

	// q1 <-- sampleFz(d, c1)
	P_sampleFz(output_2, d, c1, length);

	// compute bd^{-1}
	const std::complex<double_t> identity{1, 0};
	std::complex<double_t> * bd_1 = new std::complex<double_t>[length];
	for(uint32_t i = 0; i < length; ++i) {
		bd_1[i] = b[i]*(identity/d[i]);
	}

	// c0 := c0 + bd^{-1}(q1 - c1)
	for(uint32_t i = 0; i < length; ++i) {
		c0[i] += bd_1[i]*(output_2[i] - c1[i]);
	}

	// compute a - bd^{-1}b*
	std::complex<double_t> tmp[length];
	for(uint32_t i = 0; i < length; ++i) {
		tmp[i] = a[i] - bd_1[i]*conj(b[i]);
	}
	delete[] bd_1;

	// q0 <-- sampleFz(a - bd^{-1}b*, c0)
	P_sampleFz(output_1, tmp, c0, length);
}


void Gaussian::GaussianSampling::encapsulation::P_samplePz(Poly_t * output) const noexcept {
	// Pool of threads
	std::thread pool[THREADS_NUMBER];

	// Thread shares
	uint32_t thread_shares[THREADS_NUMBER]; 
	const uint32_t thread_share = param/THREADS_NUMBER;
	std::fill_n(thread_shares, THREADS_NUMBER, thread_share);
	for(uint32_t i = 0; i < param - thread_share*THREADS_NUMBER; ++i) {
		++thread_shares[i];
	}

	// Lambda for q
	auto lambda1 = [this](Poly_t * outputs, const uint32_t length) {
		Gauss_t noise = *(P_gaussianNoise.get());
		for(uint32_t i = 0; i < length; ++i) {
			(*outputs).set(noise);
			(*outputs).ntt_pow_phi();
			++outputs;
		}
	};

	// compute q
	Poly_t * outputIndex = output + 2;
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i] = std::thread(lambda1, outputIndex, thread_shares[i]);
		outputIndex += thread_shares[i];
	}


	// Lambda for c0 and c1
	auto lambda2 = [](Poly_t * outputs, const Poly_t * r, const Poly_t * e, const Poly_t * q, const uint32_t length) {
//		outputs[0].ntt_pow_phi();
//		outputs[1].ntt_pow_phi();
		for(uint32_t i = 0; i < length; ++i) {
			outputs[0] = outputs[0] + r[i]*q[i];
			outputs[1] = outputs[1] + e[i]*q[i]; 
		}
		outputs[0].invntt_pow_invphi();
		outputs[1].invntt_pow_invphi();
	};

	// compute c0 and c1
	Poly_t poly_c0{0};	
	Poly_t poly_c1{0};
	const Poly_t * rIndex = P_r;
	const Poly_t * eIndex = P_e;
	const Poly_t * qIndex = output + 2;
	Poly_t * poly_cs = new Poly_t[2*THREADS_NUMBER]{0};
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i].join();
	}
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		const uint32_t share = thread_shares[i];
		pool[i] = std::thread(lambda2, poly_cs + i*2, rIndex, eIndex, qIndex, share);
		rIndex += share;
		eIndex += share;
		qIndex += share;

	}
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i].join();
	}
//	poly_c0.ntt_pow_phi();
//	poly_c1.ntt_pow_phi();
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		poly_c0 = poly_c0 + *(poly_cs + i*2);
		poly_c1 = poly_c1 + *(poly_cs + 1 + i*2);
	}
	delete[] poly_cs;

	// get c0 and c1 coefs
//	poly_c0.invntt_pow_invphi();
//	poly_c1.invntt_pow_invphi();
//	auto c0Coefs = poly_c0.poly2mpz();
//	auto c1Coefs = poly_c1.poly2mpz();
//	const uint32_t coefsSize = c0Coefs.size();


	const double_t scalar = P_z2;
	std::complex<double_t> * fft_c0 = new std::complex<double_t>[dimension];
	std::complex<double_t> * fft_c1 = new std::complex<double_t>[dimension];
	auto lambda3 = [this, scalar](std::complex<double_t> * fft_c, const Poly_t & poly_c) {
		std::complex<double_t> * c = new std::complex<double_t>[dimension];
		uint32_t index = 0;
		for(auto & coef : poly_c.poly_obj()) {
			c[index++] = std::complex<double_t>{scalar*coef, 0};
		}
		fastFft.fftForward(fft_c, c);
		delete[] c;
	};
	pool[0] = std::thread(lambda3, fft_c0, poly_c0);
	pool[1] = std::thread(lambda3, fft_c1, poly_c1);
	pool[0].join();
	pool[1].join();

	/*
	for(auto & coef : poly_c0.poly_obj()) {
		c0[index++] = std::complex<double_t>{scalar*coef, 0};
//		c0[i] = std::complex<double_t>{scalar*mpz_get_ui(c0Coefs[i]), 0};
//		c1[i] = std::complex<double_t>{scalar*mpz_get_ui(c1Coefs[i]), 0};
	}
	index = 0;
	for(auto & coef : poly_c1.poly_obj()) {
		c1[index++] = std::complex<double_t>{scalar*coef, 0};
	}

	const uint32_t coefsSize = index;
	
	fastFft.fftForward(fft_c0, c0);
	fastFft.fftForward(fft_c1, c1);
	*/

	// compute c
	const uint32_t coefsSize = dimension;
	std::complex<double_t> * fft_c = new std::complex<double_t>[coefsSize << 1];
	fastFft.ffmerge(fft_c, fft_c0, fft_c1, coefsSize);
//	delete[] c0;
//	delete[] c1;
	delete[] fft_c0;
	delete[] fft_c1;


	// compute p0 and p1
	const uint32_t n = coefsSize;
//	std::complex<double_t> fft_p0[n];
//	std::complex<double_t> fft_p1[n];
	std::complex<double_t> * fft_p0 = new std::complex<double_t>[n];
	std::complex<double_t> * fft_p1 = new std::complex<double_t>[n];
	P_sample2z(fft_p0, fft_p1, P_a, P_b, P_d, fft_c, n);


	// fft inverse
	auto lambda4 = [n, this](Poly_t * outputs, std::complex<double_t> * fft_p) {
		std::complex<double_t> ifft_p[n];
		fastFft.fftInverse(ifft_p, fft_p);
//		delete[] fft_p;
		
		// get p coefs
		NFL_POLY_COEF_TYPE coefs_p[n];
		for(uint32_t i = 0; i < n; ++i) {
			coefs_p[i] = (uint32_t) ifft_p[i].real();
		}
		outputs->set(coefs_p, coefs_p + n);
		outputs->ntt_pow_phi();
	};

	pool[0] = std::thread(lambda4, output, fft_p0);
	pool[1] = std::thread(lambda4, output + 1, fft_p1);
	pool[0].join();
	pool[1].join();

	/*
	std::complex<double_t> ifft_p0[n];
	std::complex<double_t> ifft_p1[n];
	fastFft.fftInverse(ifft_p0, fft_p0);
	fastFft.fftInverse(ifft_p1, fft_p1);

	delete[] fft_p0;
	delete[] fft_p1;

	// get p0 and p1 coefs
	NFL_POLY_COEF_TYPE coefs_p0[n];
	NFL_POLY_COEF_TYPE coefs_p1[n];
	for(uint32_t i = 0; i < n; ++i) {
		coefs_p0[i] = (uint32_t) ifft_p0[i].real();
		coefs_p1[i] = (uint32_t) ifft_p1[i].real();
	}

	// build the polynom and ntt
	output[0].set(coefs_p0, coefs_p0 + n);
	output[0].ntt_pow_phi();
	output[1].set(coefs_p1, coefs_p1 + n);
	output[1].ntt_pow_phi();
	*/
}



void Gaussian::GaussianSampling::encapsulation::P_preCompute(void) const noexcept {
	// compute SamplePz
	Poly_t * pSamples{new Poly_t[param + 2]};
	P_samplePz(pSamples);

	// add to the precomputed perturbations
	std::lock_guard<std::mutex> lock(P_precomputed_mutex);
	P_precomputedPerturbations.push(pSamples);
}



/* --------------------------------------------------------------------------------------------------- */
/* Private Interface --------------------------------------------------------------------------------- */
void Gaussian::GaussianSampling::G_sampleGsMonoThread(int32_t * output, const uint64_t * coset, const uint32_t length) const noexcept {
	for(uint32_t i = 0; i < length; ++i) {
		impl->G_sampleG(output, coset[i]);
		output += impl->param;
	}
}

/* Private Interface End ----------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------------------------- */


/* --------------------------------------------------------------------------------------------------- */
/* Public Interface ---------------------------------------------------------------------------------- */
/* @Override */
Gaussian::GaussianSampling::GaussianSampling(void) noexcept {}

/* @Override */
Gaussian::GaussianSampling::GaussianSampling(const uint32_t n, const uint32_t k, const uint64_t q, const double_t s, const double_t ps, const double_t alpha, const uint32_t lambda, const Poly_t * r, const Poly_t * e) noexcept : \
impl(new encapsulation {.dimension = n, .param = k}) {
	impl->sigma = s/3;

	// get all the bits of q
	getBits(impl->modulus, q, k);

	// random generator for sampleZ
	std::random_device rd;
	impl->randomGenerator_t.reset(new std::mt19937_64(rd()));
	impl->randomGenerator = impl->randomGenerator_t.get();



	// compute G_l, G_h and G_d
	impl->G_l[0] = sqrt(2*(1 + 1/k) + 1);
	impl->G_h[0] = 0;
	for(uint32_t i = 1; i < k; ++i) {
		impl->G_l[i] = sqrt(2*(1 + 1/(k - i)));
		impl->G_h[i] = sqrt(2*(1 - 1/(k - i + 1)));
	}
	impl->G_d[0] = impl->modulus[0]/2;
	for(uint32_t i = 1; i < k; ++i) {
		impl->G_d[i] = (impl->G_d[i - 1] + impl->modulus[i])/2;
	}


	// compute P_z0, P_z1 and P_z2
	const double_t ps_2 = ps*ps;
	const double_t alpha_2 = alpha*alpha;
	impl->P_z0 = (ps_2 - alpha_2)/(alpha_2*ps_2);
	impl->P_z1 = sqrt(ps_2 - alpha_2);
	impl->P_z2 = (-alpha_2)/(ps_2 - alpha_2);


	// prebuild the Gaussian noise for P
	impl->P_fastGaussNoise.reset(new fastGauss_t(5, lambda, n));

//	impl->P_fastGaussNoice.reset(new fastGauss_t(impl->P_z1, lambda, n));
	impl->P_gaussianNoise.reset(new Gauss_t(impl->P_fastGaussNoise.get()));


	// store r and e
	impl->P_r = r;
	impl->P_e = e;	


	// compute P_a, P_b and P_d
	Poly_t rr{0};
	Poly_t re{0};
	Poly_t ee{0};
//	rr.ntt_pow_phi();
//	re.ntt_pow_phi();
//	ee.ntt_pow_phi();
	for(uint32_t i = 0; i < k; ++i) {
		rr = rr + r[i]*r[i];
		re = re + r[i]*e[i];
		ee = ee + e[i]*e[i];
	}
	rr.invntt_pow_invphi();
	re.invntt_pow_invphi();
	ee.invntt_pow_invphi();

	auto rrCoefs = rr.poly2mpz();
	auto reCoefs = re.poly2mpz();
	auto eeCoefs = ee.poly2mpz();
	const uint32_t coefsSize = rrCoefs.size();


	std::complex<double_t> tmp_a[coefsSize];
	std::complex<double_t> tmp_b[coefsSize];
	std::complex<double_t> tmp_d[coefsSize];
	for(uint32_t i = 0; i < coefsSize; ++i) {
		tmp_a[i] = std::complex<double_t>(ps_2 - impl->P_z0*mpz_get_ui(rrCoefs[i]), 0);
		tmp_b[i] = std::complex<double_t>(-impl->P_z0*mpz_get_ui(reCoefs[i]), 0);
		tmp_d[i] = std::complex<double_t>(ps_2 - impl->P_z0*mpz_get_ui(eeCoefs[i]), 0);
	}

	
	// compute the fft for P_a, P_b and P_d
	impl->fastFft.fftForward(impl->P_a, tmp_a);
	impl->fastFft.fftForward(impl->P_b, tmp_b);
	impl->fastFft.fftForward(impl->P_d, tmp_d);

	// precompute
//	preCompute(5);
}

/* @Override */
Gaussian::GaussianSampling::GaussianSampling(GaussianSampling && other) noexcept {
	impl.reset(other.impl.get());
}

/* @Override */
Gaussian::GaussianSampling::~GaussianSampling(void) noexcept {
}


/* @Override */
Gaussian::GaussianSampling & Gaussian::GaussianSampling::operator= (GaussianSampling && other) noexcept {
	impl.reset(other.impl.get());
	return *this;
}


/* @Override */
void Gaussian::GaussianSampling::preCompute(const uint32_t n) noexcept {
	// Pool of threads
	std::thread pool[THREADS_NUMBER];
	uint32_t thread_shares[THREADS_NUMBER];

	// Thread shares for G_preCompute
	const uint32_t thread_shareG = (n*impl->dimension)/THREADS_NUMBER;
	std::fill_n(thread_shares, THREADS_NUMBER, thread_shareG);
	for(uint32_t i = 0; i < n*impl->dimension - thread_shareG*THREADS_NUMBER; ++i) {
		++thread_shares[i];
	}

	// Thread function for G_preCompute
	auto lambdaG = [this](const uint32_t length) {
		for(uint32_t i = 0; i < length; ++i) {
			impl->G_preCompute();
		}
	};

	// precompute all the G perturbations
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i] = std::thread(lambdaG, thread_shares[i]);
	}
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i].join();
	}


	// P_preCompute
	const uint32_t thread_shareP = n/THREADS_NUMBER;
	std::fill_n(thread_shares, THREADS_NUMBER, thread_shareP);
	for(uint32_t i = 0; i < n - thread_shareP*THREADS_NUMBER; ++i) {
		++thread_shares[i];
	}

	// Thread function for P_preCompute
	auto lambdaP = [this](const uint32_t length) {
		for(uint32_t i = 0; i < length; ++i) {
			impl->P_preCompute();
		}
	};
	// precompute all the P perturbations
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i] = std::thread(lambdaP, thread_shares[i]);
	}
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i].join();
	}
}


/* @Override */
void Gaussian::GaussianSampling::sampleGPoly(Poly_t * polys, const uint64_t * coset) const noexcept {
	const uint32_t n = impl->dimension;
	const uint32_t k = impl->param;

	// Pool of threads
	std::thread pool[THREADS_NUMBER];

	// Thread shares
	uint32_t thread_shares[THREADS_NUMBER]; 
	const uint32_t thread_share = n/THREADS_NUMBER;
	std::fill_n(thread_shares, THREADS_NUMBER, thread_share);
	for(uint32_t i = 0; i < n - thread_share*THREADS_NUMBER; ++i) {
		++thread_shares[i];
	}

	
	// compute all the sampling
	uint32_t * samplingz = new uint32_t[n*k];
	uint32_t * samplingzInd  = samplingz;
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		const uint32_t share = thread_shares[i];
		pool[i] = std::thread(&GaussianSampling::G_sampleGsMonoThread, this, (int32_t *) samplingzInd, coset, share);
		coset += share;
		samplingzInd += (k*share);
	}
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i].join();
	}

	// transpose the samplingz matrix
	NFL_POLY_COEF_TYPE * samplingzT = new NFL_POLY_COEF_TYPE[n*k];
	const uint32_t * samplingzInd2  = samplingz;
	NFL_POLY_COEF_TYPE * samplingzTInd = samplingzT;
	for(uint32_t i = 0; i < k; ++i) {
		for(uint32_t j = 0; j < n; ++j) {
			*samplingzTInd++ = (NFL_POLY_COEF_TYPE) *(samplingzInd2 + k*j);
		}
		++samplingzInd2;
	}

	delete[] samplingz;


	// threads for polyonms and ntt
	auto lambda = [n](Poly_t * polyz, const NFL_POLY_COEF_TYPE * coefs, const uint32_t length) {
		for(uint32_t i = 0; i < length; ++i) {
			polyz[i].set(coefs, coefs + n);
			polyz[i].ntt_pow_phi();
			coefs += n;
		}
	};
	const uint32_t thread_share2 = k/THREADS_NUMBER;
	std::fill_n(thread_shares, THREADS_NUMBER, thread_share2);
	for(uint32_t i = 0; i < k - thread_share2*THREADS_NUMBER; ++i) {
		++thread_shares[i];
	}

	// compute the polynoms and ntt
	samplingzTInd = samplingzT;
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		const uint32_t share = thread_shares[i];
		pool[i] = std::thread(lambda, polys, samplingzTInd, share);
		polys += share;
		samplingzTInd += (n*share);
	}
	for(uint32_t i = 0; i < THREADS_NUMBER; ++i) {
		pool[i].join();
	}

	delete[] samplingzT;
}


/* @Override */
Poly_t * Gaussian::GaussianSampling::samplePz(void) const noexcept {
	std::lock_guard<std::mutex> lock(impl->P_precomputed_mutex);
	Poly_t * result = impl->P_precomputedPerturbations.top();
	impl->P_precomputedPerturbations.pop();
	return result;
}

/* @Override */
/*
void Gaussian::GaussianSampling::test(nfl::poly_p<NFL_POLY_COEF_TYPE, NFL_POLY_N, NFL_POLY_Q_BITS> * polys) noexcept {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(1, 65642);

    uint64_t coset[512];
    for(uint32_t i = 0; i < 512; ++i) {
    	coset[i] = (uint64_t) dis(gen);
    }

//	std::fill_n(coset, 256, 445U);
	sampleGPoly(polys, coset);
}
*/


/* Public Interface End ------------------------------------------------------------------------------ */
/* --------------------------------------------------------------------------------------------------- */


/* classical version
static int64_t sampleZ(const double_t sigma, const double_t center, const uint32_t tau = 6) {
	// compute some constants
	const double_t cst = (-1)/(2 * (sigma*sigma));
	const int64_t max = ceil(sigma*tau);

	// random engines
	std::random_device rd;
	std::mt19937_64 randomGenerator(rd());

	std::uniform_int_distribution<int64_t> uniformDist1(-max, max);
	std::uniform_real_distribution<double_t> uniformDist2(0.0, 1.0);


	// sampling by rejection
	int64_t x = 0;
	const int64_t center_z = ceil(center);
	do {
		// get a random from a uniform distribution~unif(-max, max] and center it
		x = center_z + uniformDist1(randomGenerator);

		// compute the trial
		const double_t tmp = x - center;
		const double_t z = exp(tmp*tmp*cst);

		// get a random from a uniform distribution~unif(0, 1]
		const double_t y = uniformDist2(randomGenerator);

	} while(y >= z)
	return x;
}
*/