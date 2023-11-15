#ifndef ELEM_SYM_POLY_ELEMENTARY_FFT_INCLUDED
#define ELEM_SYM_POLY_ELEMENTARY_FFT_INCLUDED

#include <stdint.h>
#include <stdio.h>
#include <common.h>

#define data_t mpfr_t
#define array_t data_t * restrict
#define matrix_t data_t ** restrict

extern void mpfr_init_look_up_table(int N, matrix_t reals, matrix_t imags);

extern void mpfr_init_look_up_inverse(int N, matrix_t reals, matrix_t imags);

extern void mpfr_free_look_up_table(int N, matrix_t reals, matrix_t imags);

extern void mpfr_recursive_fft_half_zero(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n);

extern void mpfr_recursive_inverse_fft(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n);

#undef data_t
#undef array_t
#undef matrix_t

#endif // ELEM_SYM_POLY_ELEMENTARY_FFT_INCLUDED