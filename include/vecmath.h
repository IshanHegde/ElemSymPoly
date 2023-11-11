#include <immintrin.h>


typedef struct complex_4{
    __m256d real;
    __m256d imag;
} complex_t;


static inline complex_t ADD(complex_t a, complex_t b){

    complex_t ret;
    ret.real = _mm256_add_pd(a.real, b.real);
    ret.imag = _mm256_add_pd(a.imag, b.imag);

    return ret;
} 



static inline complex_t SUB(complex_t a, complex_t b){

    complex_t ret;
    ret.real = _mm256_sub_pd(a.real, b.real);
    ret.imag = _mm256_sub_pd(a.imag, b.imag);

    return ret;
}


static inline complex_t MUL(complex_t a, complex_t b){

    complex_t ret;
    ret.real = _mm256_sub_pd(_mm256_mul_pd(a.real,b.real),_mm256_mul_pd(a.imag,b.imag));
    ret.imag = _mm256_add_pd(_mm256_mul_pd(a.real,b.imag),_mm256_mul_pd(a.imag,b.real));
    return ret;
}


static inline complex_t LOAD(double * restrict reals, double * restrict imags){

    complex_t ret;
    ret.real = _mm256_load_pd(reals);
    ret.imag = _mm256_load_pd(imags);

    return ret;
}



static inline void STORE(double * restrict reals, double * restrict imags, complex_t val){

    _mm256_store_pd(reals, val.real);
    _mm256_store_pd(imags, val.imag);

}