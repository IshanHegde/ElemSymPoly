/*
Copyright (C) 2023 Ishan Hegde

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA */

#ifdef __AVX2__
	#include <immintrin.h>

	typedef struct complex_4{
		__m256d real;
		__m256d imag;
	} complex_t;

	#define SIMD_ADD(a,b) _mm256_add_pd((a),(b))
	#define SIMD_SUB(a,b) _mm256_sub_pd((a),(b))
	#define SIMD_MUL(a,b) _mm256_mul_pd((a),(b))
	#define SIMD_LOAD(a) _mm256_load_pd((a))
	#define SIMD_STORE(a,b) _mm256_store_pd((a),(b))

	#define LOOP_INCRIEMENT 4
#elif defined(__ARM_NEON)
	#include <arm_neon.h>

	typedef struct complex_4 {
		float64x2_t real;
		float64x2_t imag;
	} complex_t;

	#define SIMD_ADD(a,b) vaddq_f64((a),(b))
	#define SIMD_SUB(a,b) vsubq_f64((a),(b))
	#define SIMD_MUL(a,b) vmulq_f64((a),(b))
	#define SIMD_LOAD(a) vld1q_f64((a))
	#define SIMD_STORE(a,b) vst1q_f64((a),(b))

	#define LOOP_INCRIEMENT 2
#else
	typedef struct complex_4 {
		double real;
		double imag;
	} complex_t;
	#define SIMD_ADD(a,b) (a+b)
	#define SIMD_SUB(a,b) (a-b)
	#define SIMD_MUL(a,b) (a*b)
	#define SIMD_LOAD(a) (a)
	#define SIMD_STORE(a,b) (a=b)

	#define LOOP_INCRIEMENT 1
#endif

static inline complex_t ADD(complex_t a, complex_t b){

    complex_t ret;
    ret.real = SIMD_ADD(a.real, b.real);
    ret.imag = SIMD_ADD(a.imag, b.imag);

    return ret;
} 



static inline complex_t SUB(complex_t a, complex_t b){

    complex_t ret;
    ret.real = SIMD_SUB(a.real, b.real);
    ret.imag = SIMD_SUB(a.imag, b.imag);

    return ret;
}


static inline complex_t MUL(complex_t a, complex_t b){

    complex_t ret;
    ret.real = SIMD_SUB(SIMD_MUL(a.real,b.real),SIMD_MUL(a.imag,b.imag));
    ret.imag = SIMD_ADD(SIMD_MUL(a.real,b.imag),SIMD_MUL(a.imag,b.real));
    return ret;
}


static inline complex_t LOAD(double * restrict reals, double * restrict imags){

    complex_t ret;
    ret.real = SIMD_LOAD(reals);
    ret.imag = SIMD_LOAD(imags);

    return ret;
}



static inline void STORE(double * restrict reals, double * restrict imags, complex_t val){

    SIMD_STORE(reals, val.real);
	SIMD_STORE(imags, val.imag);
}