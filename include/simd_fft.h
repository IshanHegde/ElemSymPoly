//
// Created by ishan on 11/12/23.
//

#ifndef ELEM_SYM_POLY_SIMD_FFT_H
#define ELEM_SYM_POLY_SIMD_FFT_H

#include <common.h>

#define data_t double
#define array_t data_t * restrict
#define matrix_t data_t ** restrict


extern void init_look_up_table(int N, matrix_t reals, matrix_t imags);

extern void init_look_up_inverse(int N, matrix_t reals, matrix_t imags);

extern void free_look_up_table(int N, matrix_t reals, matrix_t imags);

extern void recursive_fft_half_zero(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n);

extern void recursive_inverse_fft(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n);


#endif //ELEM_SYM_POLY_SIMD_FFT_H
