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

#include <polynomial.h>
#include <mpfr_fft.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define data_t mpfr_t
#define array_t data_t *
#define matrix_t data_t **

#define state_t poly_mul_state_t

struct state_t{

    int N;

    array_t A_reals;
	array_t A_imags;

    array_t B_reals;
	array_t B_imags;

    array_t A_out_reals;
    array_t A_out_imags;

	array_t B_out_reals;
	array_t B_out_imags;

    array_t C_out_reals;
    array_t C_out_imags;

    matrix_t w_reals;
    matrix_t w_imags;
    matrix_t w_reals_inverse;
    matrix_t w_imags_inverse;

    array_t temp_C_reals;
    array_t temp_C_imags;
};


state_t init_polynomial_mul_state(int poly_size){

    int N = 2* poly_size;

    state_t state = malloc(sizeof(struct state_t));

	state->A_reals = ALLOC(sizeof(data_t) * N);
	state->A_imags = ALLOC(sizeof(data_t) * N);

	state->B_reals = ALLOC(sizeof(data_t) * N);
	state->B_imags = ALLOC(sizeof(data_t) * N);

	state->A_out_reals = ALLOC(sizeof(data_t) * N);
	state->A_out_imags = ALLOC(sizeof(data_t) * N);

	state->B_out_reals = ALLOC(sizeof(data_t) * N);
	state->B_out_imags = ALLOC(sizeof(data_t) * N);

	state->C_out_reals = ALLOC(sizeof(data_t) * N);
	state->C_out_imags = ALLOC(sizeof(data_t) * N);

	state->temp_C_reals = ALLOC(sizeof(data_t) * N);
	state->temp_C_imags = ALLOC(sizeof(data_t) * N);

    for (int i = 0; i< N;i++){
		mpfr_inits2(PRECISION,
					state->A_reals[i],
					state->A_imags[i],
					state->B_reals[i],
					state->B_imags[i],
					state->A_out_reals[i],
					state->A_out_imags[i],
					state->B_out_reals[i],
					state->B_out_imags[i],
					state->C_out_reals[i],
					state->C_out_imags[i],
					state->temp_C_reals[i],
					state->temp_C_imags[i],
					(mpfr_ptr)NULL);

		mpfr_set_d(state->A_imags[i], 0, MPFR_RNDN);
		mpfr_set_d(state->B_imags[i], 0, MPFR_RNDN);
    }
    //memset(state->C_out_reals, 0, sizeof(data_t) * N);
    //memset(state->C_out_imags, 0,  sizeof(data_t) * N);

    state->w_reals = ALLOC(sizeof(array_t) * log2(N));
    state->w_imags = ALLOC(sizeof(array_t) * log2(N));

    state->w_reals_inverse = ALLOC(sizeof(array_t) * log2(N));
    state->w_imags_inverse = ALLOC(sizeof(array_t) * log2(N));

    mpfr_init_look_up_table(N,state->w_reals,state->w_imags);
    mpfr_init_look_up_inverse(N,state->w_reals_inverse,state->w_imags_inverse);
    state->N = N;

    return state;
}


void update_polynomial_mul_state(state_t state, array_t A, array_t B, int poly_size){
    
    int N = 2* poly_size;

    state->N = N;


	for (int i = 0; i < poly_size; i++){
		mpfr_set(state->A_reals[i], A[i], MPFR_RNDN);
		mpfr_set(state->B_reals[i], B[i], MPFR_RNDN);
	}

	for (int i = poly_size; i < N; i++){
		mpfr_set_d(state->A_reals[i], 0, MPFR_RNDN);
		mpfr_set_d(state->B_reals[i], 0, MPFR_RNDN);
	}


}

array_t polynomial_multiply(state_t state){

    int N = state->N;

    mpfr_recursive_fft_half_zero(state->A_reals,state->A_imags,state->A_out_reals,state->A_out_imags, state->w_reals, state->w_imags,1,N);
    mpfr_recursive_fft_half_zero(state->B_reals,state->B_imags,state->B_out_reals,state->B_out_imags, state->w_reals, state->w_imags,1,N);


	data_t aux_0, aux_1, aux_2, aux_3;
	mpfr_inits2(PRECISION, aux_0, aux_1, aux_2, aux_3, (mpfr_ptr)NULL);
	for (int i =0; i< N;i++){
		mpfr_mul(aux_0, state->A_out_reals[i], state->B_out_reals[i], MPFR_RNDN);
		mpfr_mul(aux_1, state->A_out_imags[i], state->B_out_imags[i], MPFR_RNDN);
		mpfr_sub(state->temp_C_reals[i],aux_0,aux_1, MPFR_RNDN);

		mpfr_mul(aux_2, state->A_out_reals[i], state->B_out_imags[i], MPFR_RNDN);
		mpfr_mul(aux_3, state->A_out_imags[i], state->B_out_reals[i], MPFR_RNDN);

		mpfr_add(state->temp_C_imags[i],aux_2,aux_3, MPFR_RNDN);
	}
	mpfr_clears(aux_0, aux_1, aux_2, aux_3, (mpfr_ptr)NULL);

    mpfr_recursive_inverse_fft(state->temp_C_reals,state->temp_C_imags,state->C_out_reals,state->C_out_imags,state->w_reals_inverse,state->w_imags_inverse,1,N);

    for (int i = 0;i < N;i++){
        mpfr_div_ui(state->C_out_reals[i],state->C_out_reals[i],N,MPFR_RNDN);
    }


    return state->C_out_reals;
}

void free_polynomial_mul_state(state_t state){

    for (int i =0;i < state->N;i++){
        mpfr_clears(state->A_reals[i],
                    state->A_imags[i],
                    state->B_reals[i],
                    state->B_imags[i],
                    state->A_out_reals[i],
                    state->A_out_imags[i],
                    state->B_out_reals[i],
                    state->B_out_imags[i],
                    state->temp_C_reals[i],
                    state->temp_C_imags[i],
                    state->C_out_reals[i],
                    state->C_out_imags[i],
                    (mpfr_ptr)NULL);
    }
    free(state->A_reals);
    free(state->A_imags);
    free(state->B_reals);
    free(state->B_imags);

    free(state->A_out_reals);
    free(state->A_out_imags);
    free(state->B_out_reals);
    free(state->B_out_imags);

    free(state->temp_C_reals);
    free(state->temp_C_imags);
    free(state->C_out_reals);
    free(state->C_out_imags);


    mpfr_free_look_up_table(state->N,state->w_reals,state->w_imags);
    mpfr_free_look_up_table(state->N,state->w_reals_inverse,state->w_imags_inverse);

    free(state);
}

#undef state_t
#undef data_t
#undef array_t
#undef matrix_t