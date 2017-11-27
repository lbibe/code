#ifndef FAST_FFT_HPP
#define FAST_FFT_HPP

#include<cstdint>
#include<cmath>
#include<complex>
#include<utility>


namespace FFT {
	template<unsigned N, typename T = double_t>
	class DanielsonLanczos {
		private:
			const DanielsonLanczos<N/2, T> next;

		public:
			void apply(T * data) const noexcept {
				next.apply(data);
				next.apply(data + N);

				T wtemp = sin(M_PI/N);
				T wpr = -2.0*wtemp*wtemp;
				T wpi = -sin(2*M_PI/N);
				T wr = 1.0;
				T wi = 0.0;

				T tempr, tempi;
				for (uint32_t i = 0; i < N; i+=2) {
					tempr = data[i+N]*wr - data[i+N+1]*wi;
					tempi = data[i+N]*wi + data[i+N+1]*wr;
					data[i+N] = data[i] - tempr;
					data[i+N+1] = data[i+1] - tempi;
					data[i] += tempr;
					data[i+1] += tempi;

					wtemp = wr;
					wr += wr*wpr - wi*wpi;
					wi += wi*wpr + wtemp*wpi;
				}
			}
	};


	template<typename T>
	class DanielsonLanczos<1, T> {
		public:
	   		void apply(T * data) const noexcept {}
	};


	template<unsigned N, typename T = double_t>
	class FastFFT {
		public:
			/** 
			 * @Constructor
			 */
			FastFFT(void) noexcept : recursion{} {
				const double_t theta    = -2.*M_PI/(2*N);
				const double_t invTheta =  2.*M_PI/(2*N);

				for(uint32_t i = 0; i < N; ++i) {
					twiddles[i]    = exp(std::complex<T>(0, theta*i));
					invTwiddles[i] = exp(std::complex<T>(0, invTheta*i));
				}
			}

			/** 
			 * @Destructor
			 */
			~FastFFT(void) noexcept {}



			/** 
			 * @
			 */
			void fftForward(std::complex<T> * output, const std::complex<T> * input) const noexcept {
				// prepare data for fft
				T data[N<<1];
				for(uint32_t i = 0, j = 0; j < N; ++j) {
					data[i++] = input[j].real();
					data[i++] = input[j].imag();
				}
				fft(data);

				// prepare output data
				for(uint32_t i = 0, j = 0; j < N; ++j, i += 2) {
					output[j] = std::complex<T>{data[i], data[i + 1]};
				}
			}
			


			/** 
			 * @
			 */
			void fftInverse(std::complex<T> * output, const std::complex<T> * input) const noexcept {
				// conjugate
				for(uint32_t i = 0; i < N; ++i) {
					output[i] = std::conj(input[i]);
				}

				std::complex<T> intermediate[N];
				fftForward(intermediate, output);

				// conjugate and scale
				const T scalar = ((T) 1)/N;
				for(uint32_t i = 0; i < N; ++i) {
					output[i] = scalar*std::conj(intermediate[i]);
				}

			}

			/**
			 * @
			 */
			void ffmerge(std::complex<T> * output, const std::complex<T> * input_1, const std::complex<T> * input_2, const uint32_t length) const noexcept {
				const uint32_t multiplicator = N/length;
				for(uint32_t i = 0; i < length; ++i) {
					const std::complex<T> e = input_1[i];
            		const std::complex<T> o = input_2[i];
					output[i]          = e + o*twiddles[multiplicator*i];
					output[i + length] = e - o*twiddles[multiplicator*i];
				}
			}


			/**
			 * @
			 */
			void ffsplit(std::complex<T> * output_1, std::complex<T> * output_2, const std::complex<T> * input, const uint32_t length) const noexcept {
				const uint32_t d             = length >> 1;
				const uint32_t multiplicator = N/d;
				for(uint32_t i = 0; i < d; ++i) {
					const std::complex<T> e = input[i];
            		const std::complex<T> o = input[i + d];
					output_1[i] = 0.5*(e + o);
					output_2[i] = 0.5*(e - o)*invTwiddles[multiplicator*i];
				}
			}


		private:
			std::complex<T> twiddles[N];
			std::complex<T> invTwiddles[N];
			const DanielsonLanczos<N, T> recursion;

			void scramble(double_t * data, const uint32_t length) const noexcept {
				uint32_t j = 1;
				for(uint32_t i = 1; i < length<<1; i+=2) {
					if(j > i) {
						std::swap(data[j-1], data[i-1]);
						std::swap(data[j], data[i]);
					}
					uint32_t m = length;
					while((m >= 2) && (j > m)) {
						j -= m;
						m >>= 1;
					}
					j += m;
				}
			}

			void fft(T * data) const noexcept {
				scramble(data, N);
				recursion.apply(data);
			}
	};
}









#endif