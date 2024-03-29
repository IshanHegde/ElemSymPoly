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

#ifndef ELEM_SYM_POLY_ELEMENTARY_FFT_INCLUDED
#define ELEM_SYM_POLY_ELEMENTARY_FFT_INCLUDED

#include <stdint.h>
#include <stdio.h>
#include "common.h"

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