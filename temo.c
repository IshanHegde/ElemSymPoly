#include <immintrin.h>

typedef struct complex_4{
    __m256d real;
    __m256d imag;
} complex_4;

static inline complex_4 ADD_4(complex_4 a, complex_4 b){

    complex_4 ret;
    ret.real = _mm256_add_pd(a.real, b.real);
    ret.imag = _mm256_add_pd(a.imag, b.imag);

    return ret;
} 

static inline complex_4 LOAD_4(double * restrict reals, double * restrict imags){

    complex_4 ret;
    ret.real = _mm256_load_pd(reals);
    ret.imag = _mm256_load_pd(imags);

    return ret;
}
static inline void STORE_4(double * restrict reals, double * restrict imags, complex_4 val){

    _mm256_store_pd(reals, val.real);
    _mm256_store_pd(imags, val.imag);

}

static inline complex_4 MUL_4(complex_4 a, complex_4 b){

    complex_4 ret;
    ret.real = _mm256_sub_pd(_mm256_mul_pd(a.real,b.real),_mm256_mul_pd(a.imag,b.imag));
    ret.imag = _mm256_add_pd(_mm256_mul_pd(a.real,b.imag),_mm256_mul_pd(a.imag,b.real));

    return ret;
}

static inline complex_4 SUB_4(complex_4 a, complex_4 b){

    complex_4 ret;
    ret.real = _mm256_sub_pd(a.real, b.real);
    ret.imag = _mm256_sub_pd(a.imag, b.imag);

    return ret;
}


