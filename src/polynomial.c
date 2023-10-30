#include <polynomial.h>
#include <dft.h>
#include <string.h>

#define data_t double
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
    int alignment = 32;

    state_t state = malloc(sizeof(struct state_t));
    
    ALLOC_ALIGNED(state->A_reals, alignment, sizeof(data_t) * N);

    ALLOC_ALIGNED(state->B_reals, alignment, sizeof(data_t) * N);

    ALLOC_ALIGNED(state->temp_C_reals, alignment, sizeof(data_t) * N);
    ALLOC_ALIGNED(state->temp_C_imags, alignment, sizeof(data_t) * N);
    
    ALLOC_ALIGNED(state->A_out_reals, alignment, sizeof(data_t) * N);
    ALLOC_ALIGNED(state->A_out_imags, alignment, sizeof(data_t) * N);


    ALLOC_ALIGNED(state->C_out_reals, alignment, sizeof(data_t) * N);
    ALLOC_ALIGNED(state->C_out_imags, alignment, sizeof(data_t) * N);

    memset(state->C_out_reals, 0, sizeof(data_t) * N);
    memset(state->C_out_imags, 0,  sizeof(data_t) * N);

    state->w_reals = (matrix_t) malloc(sizeof(array_t) * log2(N));
    state->w_imags = (matrix_t) malloc(sizeof(array_t) * log2(N));

    state->w_reals_inverse = (matrix_t) malloc(sizeof(array_t) * log2(N));
    state->w_imags_inverse = (matrix_t) malloc(sizeof(array_t) * log2(N));

    init_look_up_table_d(N,state->w_reals,state->w_imags);
    init_look_up_inverse_d(N,state->w_reals_inverse,state->w_imags_inverse);
    state->N = N;

    return state;
}


void update_polynomial_mul_state(state_t state, array_t A, array_t B, int poly_size){
    
    int N = 2* poly_size;

    state->N = N;

    size_t byte_size = poly_size * sizeof(data_t);

    memcpy(state->A_reals, A, byte_size);
    memcpy(state->B_reals, B, byte_size);

    //memset(state->A_out_reals, 0, 2 * byte_size);
    //memset(state->B_out_reals, 0, 2 * byte_size);

    memset(state->A_reals + poly_size, 0, byte_size);
    memset(state->B_reals + poly_size, 0, byte_size);

    //memset(state->C_out_reals, 0, 2* byte_size);
    //memset(state->C_out_imags, 0, 2* byte_size);
}

array_t polynomial_multiply(state_t state){

    int N = state->N;

    recursive_rfft_half_zero_d(state->A_reals,state->B_reals,state->A_out_reals,state->A_out_imags,state->w_reals,state->w_imags,1,N);

    state->temp_C_reals[0] = state->A_out_reals[0] * state->A_out_imags[0];
    state->temp_C_imags[0] = 0;

    for (int i = 1;i < N;i++){
        
        data_t a = state->A_out_reals[i];
        data_t b = state->A_out_imags[i];

        data_t c = state->A_out_reals[N - i];
        data_t d = state->A_out_imags[N - i];

        state->temp_C_reals[i] = 0.5 * (a * b +  c * d);
        state->temp_C_imags[i] = 0.25 * (c * c + b * b - a * a - d * d);
    }


    recursive_inverse_fft_d(state->temp_C_reals,state->temp_C_imags,state->C_out_reals,state->C_out_imags,state->w_reals_inverse,state->w_imags_inverse,1,N);

    for (int i = 0;i < N;i++){
        state->C_out_reals[i] /= N;
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
    free(state->w_reals);
    free(state->w_imags);
    free(state->w_reals_inverse);
    free(state->w_imags_inverse);
    free(state);
}

#undef state_t
#undef data_t
#undef array_t
#undef matrix_t