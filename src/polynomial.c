#include <polynomial.h>
#include <fft.h>
#include <string.h>

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
    }
    //memset(state->C_out_reals, 0, sizeof(data_t) * N);
    //memset(state->C_out_imags, 0,  sizeof(data_t) * N);

    state->w_reals = ALLOC(sizeof(array_t) * log2(N));
    state->w_imags = ALLOC(sizeof(array_t) * log2(N));

    state->w_reals_inverse = ALLOC(sizeof(array_t) * log2(N));
    state->w_imags_inverse = ALLOC(sizeof(array_t) * log2(N));

    init_look_up_table(N,state->w_reals,state->w_imags);
    init_look_up_inverse(N,state->w_reals_inverse,state->w_imags_inverse);
    state->N = N;

    return state;
}


void update_polynomial_mul_state(state_t state, array_t A, array_t B, int poly_size){
    
    int N = 2* poly_size;

    state->N = N;

    //size_t byte_size = poly_size * sizeof(data_t);

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
	puts("input A, B:");
	for (int i =0;i < N;i++){
		mpfr_printf("%Rf +I %Rf\n ", state->A_reals[i], state->B_reals[i]);
	}
	recursive_rfft_half_zero_safe(state->A_reals,state->A_out_imags,state->A_out_reals,state->A_out_imags, state->w_reals, state->w_imags,1,N);
	recursive_rfft_half_zero_safe(state->B_reals,state->B_out_imags,state->B_out_reals,state->B_out_imags, state->w_reals, state->w_imags,1,N);


	puts("A_out:");
	for (int i =0;i < N;i++){
		mpfr_printf("%Rf +I %Rf\n ", state->A_out_reals[i], state->A_out_imags[i]);

	}

    //state->temp_C_reals[0] = state->A_out_reals[0] * state->A_out_imags[0];
	mpfr_mul(state->temp_C_reals[0], state->A_out_reals[0], state->A_out_imags[0], MPFR_RNDN);

	mpfr_t half, quarter, aux_0, aux_1;
	mpfr_init2(half, PRECISION);
	mpfr_init2(quarter, PRECISION);
    mpfr_init2(aux_0,PRECISION);
    mpfr_init2(aux_1,PRECISION);
	mpfr_set_d(half, 0.5, MPFR_RNDN);
	mpfr_set_d(quarter, 0.25, MPFR_RNDN);
    for (int i = 1;i < N;i++){
        /*
        data_t a = state->A_out_reals[i];
        data_t b = state->A_out_imags[i];

        data_t c = state->A_out_reals[N - i];
        data_t d = state->A_out_imags[N - i];
		*/
        mpfr_fmma(aux_0,state->A_out_reals[i],state->A_out_imags[i], state->A_out_reals[N-i],state->A_out_imags[N-i],MPFR_RNDN);
        mpfr_mul(state->temp_C_reals[i],half,aux_0,MPFR_RNDN);
        //state->temp_C_reals[i] = 0.5 * (a * b +  c * d);
        mpfr_fmma(aux_0,state->A_out_reals[N-i],state->A_out_reals[N-i],state->A_out_imags[i],state->A_out_imags[i],MPFR_RNDN);
        mpfr_fmma(aux_1,state->A_out_reals[i],state->A_out_reals[i],state->A_out_imags[N-i],state->A_out_imags[N-i],MPFR_RNDN);
        mpfr_add(state->temp_C_imags[i],aux_1, aux_0,MPFR_RNDN);
        mpfr_mul(state->temp_C_imags[i],state->temp_C_imags[i],quarter,MPFR_RNDN);
        //state->temp_C_imags[i] = 0.25 * (c * c + b * b - a * a - d * d);

    }
    mpfr_clear(half);
    mpfr_clear(quarter);
    mpfr_clear(aux_0);
    mpfr_clear(aux_1);

	recursive_inverse_fft_safe(state->temp_C_reals,state->temp_C_imags,state->C_out_reals,state->C_out_imags,state->w_reals_inverse,state->w_imags_inverse,1,N);

    for (int i = 0;i < N;i++){
        mpfr_div_ui(state->C_out_reals[i],state->C_out_reals[i],N,MPFR_RNDN);
    }


    return state->C_out_reals;
}

void free_polynomial_mul_state(state_t state){

    free(state->A_reals);
    free(state->B_reals);
    free(state->A_out_reals);
    free(state->A_out_imags);
    free(state->temp_C_reals);
    free(state->temp_C_imags);
    free(state->C_out_reals);
    free(state->C_out_imags);

    free_look_up_table(state->N,state->w_reals,state->w_imags);
    free_look_up_table(state->N,state->w_reals_inverse,state->w_imags_inverse);

    free(state);
}

#undef state_t
#undef data_t
#undef array_t
#undef matrix_t