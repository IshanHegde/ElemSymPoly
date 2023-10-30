#ifndef PYRASCH_FFT_INCLUDED
#define PYRASCH_FFT_INCLUDED

#include <stdint.h>
#include <stdio.h>


#define data_t double
#define array_t data_t * restrict
#define matrix_t data_t ** restrict

void init_look_up_table_d(int N, matrix_t reals, matrix_t imags);

void init_look_up_inverse_d(int N, matrix_t reals, matrix_t imags);

void recursive_fft_d(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n);

void recursive_inverse_fft_d(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n);

void recursive_rfft_half_zero_d(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n);

#undef data_t
#undef array_t
#undef matrix_t

#endif // PYRASCH_FFT_INCLUDED