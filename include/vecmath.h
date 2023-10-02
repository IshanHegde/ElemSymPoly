#include <immintrin.h>

typedef struct complex_8{

    __m256 real;
    __m256 imag;
} complex_8;


static inline complex_8 ADD(complex_8 a, complex_8 b){

    complex_8 ret;
    ret.real = _mm256_add_ps(a.real, b.real);
    ret.imag = _mm256_add_ps(a.imag, b.imag);

    return ret;
}

static inline complex_8 SUB(complex_8 a, complex_8 b){

    complex_8 ret;
    ret.real = _mm256_sub_ps(a.real, b.real);
    ret.imag = _mm256_sub_ps(a.imag, b.imag);

    return ret;
}

static inline complex_8 MUL(complex_8 a, complex_8 b){

    complex_8 ret;
    ret.real = _mm256_sub_ps(_mm256_mul_ps(a.real,b.real),_mm256_mul_ps(a.imag,b.imag));
    ret.imag = _mm256_add_ps(_mm256_mul_ps(a.real,b.imag),_mm256_mul_ps(a.imag,b.real));

    return ret;
}

static inline complex_8 LOAD(float * restrict reals, float * restrict imags){

    complex_8 ret;
    ret.real = _mm256_loadu_ps(reals);
    ret.imag = _mm256_loadu_ps(imags);

    return ret;
}

static inline void STORE(float * restrict reals, float * restrict imags, complex_8 val){

    _mm256_storeu_ps(reals, val.real);
    _mm256_storeu_ps(imags, val.imag);

}