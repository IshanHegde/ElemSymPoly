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

#include <mpfr_fft.h>
#include <math.h>
#include <stdlib.h>
#include <common.h>


#define data_t mpfr_t
#define array_t data_t * restrict
#define matrix_t data_t ** restrict



void mpfr_init_look_up_table(int N, matrix_t reals, matrix_t imags){

	int outer_array_size = (int)log2(N);

	for (int i =0; i < outer_array_size;i++){

		int inner_array_size = (int) pow(2,i);
		int aux_size = inner_array_size * 2;
		reals[i] = ALLOC(sizeof(data_t) * inner_array_size);
		imags[i] = ALLOC(sizeof(data_t) * inner_array_size);

		data_t op, aux_0, aux_1, pi;
		mpfr_inits2(PRECISION, op, aux_0, aux_1, pi, (mpfr_ptr)NULL );

		mpfr_const_pi(pi, MPFR_RNDN);
		mpfr_div_ui(aux_0, pi, aux_size, MPFR_RNDN);
		mpfr_mul_ui(aux_1, aux_0, 2, MPFR_RNDN);

		for (int j =0;j < inner_array_size;j++){

			mpfr_init2(reals[i][j], PRECISION);
			mpfr_init2(imags[i][j], PRECISION);
			mpfr_mul_ui(op, aux_1, j, MPFR_RNDN);

			mpfr_cos(reals[i][j],op, MPFR_RNDN);
			mpfr_sin(imags[i][j],op, MPFR_RNDN);

		}

		mpfr_clears(op, aux_0, aux_1, pi, (mpfr_ptr)NULL);
	}
}

void mpfr_init_look_up_inverse(int N, matrix_t reals, matrix_t imags){

	int outer_array_size = (int)log2(N);

	for (int i =0; i < outer_array_size;i++){

		int inner_array_size = (int) pow(2,i);
		int aux_size = inner_array_size * 2;
		reals[i] = ALLOC(sizeof(data_t) * inner_array_size);
		imags[i] = ALLOC(sizeof(data_t) * inner_array_size);

		data_t op, aux_0, aux_1, pi;
		mpfr_inits2(PRECISION, op, aux_0, aux_1, pi, (mpfr_ptr)NULL );

		mpfr_const_pi(pi, MPFR_RNDN);
		mpfr_div_ui(aux_0, pi, aux_size, MPFR_RNDN);
		mpfr_mul_si(aux_1, aux_0, -2, MPFR_RNDN);

		for (int j =0;j < inner_array_size;j++){

			mpfr_init2(reals[i][j], PRECISION);
			mpfr_init2(imags[i][j], PRECISION);
			mpfr_mul_ui(op, aux_1, j, MPFR_RNDN);

			mpfr_cos(reals[i][j],op, MPFR_RNDN);
			mpfr_sin(imags[i][j],op, MPFR_RNDN);

		}

		mpfr_clears(op, aux_0, aux_1, pi, (mpfr_ptr)NULL);
	}

}

void mpfr_free_look_up_table(int N, matrix_t reals, matrix_t imags){

    int outer_array_size = (int)log2(N);

    for (int i = 0; i < outer_array_size; i++){
        int inner_array_size = (int) pow(2,i);
        for (int j = 0; j < inner_array_size; j++){
            mpfr_clears(reals[i][j], imags[i][j], (mpfr_ptr)NULL);
        }
        free(reals[i]);
        free(imags[i]);
    }

    free(reals);
    free(imags);
}

void mpfr_recursive_fft_half_zero(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n){

	if (n == 4){

		mpfr_add(out_reals[0],in_reals[0],in_reals[stride], MPFR_RNDN);
		mpfr_add(out_imags[0],in_imags[0],in_imags[stride], MPFR_RNDN);

		mpfr_sub(out_reals[1],in_reals[0],in_imags[stride], MPFR_RNDN);
		mpfr_add(out_imags[1],in_imags[0],in_reals[stride], MPFR_RNDN);

		mpfr_sub(out_reals[2],in_reals[0],in_reals[stride], MPFR_RNDN);
		mpfr_sub(out_imags[2],in_imags[0],in_imags[stride], MPFR_RNDN);

		mpfr_add(out_reals[3],in_reals[0],in_imags[stride], MPFR_RNDN);
		mpfr_sub(out_imags[3],in_imags[0],in_reals[stride], MPFR_RNDN);

	}else{
		mpfr_recursive_fft_half_zero(in_reals, in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
		mpfr_recursive_fft_half_zero(in_reals + stride, in_imags + stride, out_reals + n/2, out_imags + n/2, w_reals, w_imags, stride << 1, n >> 1);

		// twiddle factor outer array index
		int aux_num = log2(n)-1;
		mpfr_t t_real, t_imag, aux_0, aux_1, aux_2, aux_3;
		mpfr_inits2(PRECISION, t_real, t_imag, aux_0, aux_1, aux_2, aux_3, (mpfr_ptr)NULL);

		for (int k =0; k < n/2; k++){

			mpfr_set(aux_2,out_reals[k],MPFR_RNDN);
			mpfr_set(aux_3,out_imags[k],MPFR_RNDN);

			mpfr_mul(aux_0, w_reals[aux_num][k], out_reals[k + n/2], MPFR_RNDN);
			mpfr_mul(aux_1, w_imags[aux_num][k], out_imags[k + n/2], MPFR_RNDN);
			mpfr_sub(t_real, aux_0, aux_1, MPFR_RNDN);

			mpfr_mul(aux_0, w_reals[aux_num][k], out_imags[k + n/2], MPFR_RNDN);
			mpfr_mul(aux_1, w_imags[aux_num][k], out_reals[k + n/2], MPFR_RNDN);
			mpfr_add(t_imag, aux_0, aux_1, MPFR_RNDN);


			mpfr_add(out_reals[k], out_reals[k], t_real, MPFR_RNDN);
			mpfr_add(out_imags[k], out_imags[k], t_imag, MPFR_RNDN);

			mpfr_sub(out_reals[k + n/2], aux_2, t_real, MPFR_RNDN);
			mpfr_sub(out_imags[k + n/2], aux_3, t_imag, MPFR_RNDN);

		}
		mpfr_clears(t_real, t_imag, aux_0, aux_1, aux_2, aux_3, (mpfr_ptr)NULL);

	}

}

void mpfr_recursive_inverse_fft(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n){

	if (n == 4){

		mpfr_add(out_reals[0], in_reals[0], in_reals[2*stride], MPFR_RNDN);
		mpfr_add(out_reals[0], out_reals[0], in_reals[stride], MPFR_RNDN);
		mpfr_add(out_reals[0], out_reals[0], in_reals[3*stride], MPFR_RNDN);

		mpfr_add(out_imags[0], in_imags[0], in_imags[2*stride], MPFR_RNDN);
		mpfr_add(out_imags[0], out_imags[0], in_imags[stride], MPFR_RNDN);
		mpfr_add(out_imags[0], out_imags[0], in_imags[3*stride], MPFR_RNDN);

		mpfr_sub(out_reals[1], in_reals[0], in_reals[2*stride], MPFR_RNDN);
		mpfr_add(out_reals[1], out_reals[1], in_imags[stride], MPFR_RNDN);
		mpfr_sub(out_reals[1], out_reals[1], in_imags[3*stride], MPFR_RNDN);

		mpfr_sub(out_imags[1], in_imags[0], in_imags[2*stride], MPFR_RNDN);
		mpfr_sub(out_imags[1], out_imags[1], in_reals[stride], MPFR_RNDN);
		mpfr_add(out_imags[1], out_imags[1], in_reals[3*stride], MPFR_RNDN);

		mpfr_add(out_reals[2], in_reals[0], in_reals[2*stride], MPFR_RNDN);
		mpfr_sub(out_reals[2], out_reals[2], in_reals[stride], MPFR_RNDN);
		mpfr_sub(out_reals[2], out_reals[2], in_reals[3*stride], MPFR_RNDN);

		mpfr_add(out_imags[2], in_imags[0], in_imags[2*stride], MPFR_RNDN);
		mpfr_sub(out_imags[2], out_imags[2], in_imags[stride], MPFR_RNDN);
		mpfr_sub(out_imags[2], out_imags[2], in_imags[3*stride], MPFR_RNDN);

		mpfr_sub(out_reals[3], in_reals[0], in_reals[2*stride], MPFR_RNDN);
		mpfr_sub(out_reals[3], out_reals[3], in_imags[stride], MPFR_RNDN);
		mpfr_add(out_reals[3], out_reals[3], in_imags[3*stride], MPFR_RNDN);

		mpfr_sub(out_imags[3], in_imags[0], in_imags[2*stride], MPFR_RNDN);
		mpfr_add(out_imags[3], out_imags[3], in_reals[stride], MPFR_RNDN);
		mpfr_sub(out_imags[3], out_imags[3], in_reals[3*stride], MPFR_RNDN);


	}else{
		mpfr_recursive_inverse_fft(in_reals, in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
		mpfr_recursive_inverse_fft(in_reals + stride, in_imags + stride, out_reals + n/2, out_imags + n/2, w_reals, w_imags, stride << 1, n >> 1);

		int aux_num = log2(n)-1;
		mpfr_t t_real, t_imag, aux_0, aux_1, aux_2, aux_3;
		mpfr_inits2(PRECISION, t_real, t_imag, aux_0, aux_1, aux_2, aux_3, (mpfr_ptr)NULL);

		for (int k =0; k < n/2; k++){

			mpfr_set(aux_2,out_reals[k],MPFR_RNDN);
			mpfr_set(aux_3,out_imags[k],MPFR_RNDN);

			mpfr_mul(aux_0, w_reals[aux_num][k], out_reals[k + n/2], MPFR_RNDN);
			mpfr_mul(aux_1, w_imags[aux_num][k], out_imags[k + n/2], MPFR_RNDN);
			mpfr_sub(t_real, aux_0, aux_1, MPFR_RNDN);

			mpfr_mul(aux_0, w_reals[aux_num][k], out_imags[k + n/2], MPFR_RNDN);
			mpfr_mul(aux_1, w_imags[aux_num][k], out_reals[k + n/2], MPFR_RNDN);
			mpfr_add(t_imag, aux_0, aux_1, MPFR_RNDN);


			mpfr_add(out_reals[k], aux_2, t_real, MPFR_RNDN);
			mpfr_add(out_imags[k], aux_3, t_imag, MPFR_RNDN);

			mpfr_sub(out_reals[k + n/2], aux_2, t_real, MPFR_RNDN);
			mpfr_sub(out_imags[k + n/2], aux_3, t_imag, MPFR_RNDN);


		}
		mpfr_clears( t_real, t_imag, aux_0, aux_1, aux_2, aux_3, (mpfr_ptr)NULL);

	}
}

#undef data_t
#undef array_t
#undef matrix_t