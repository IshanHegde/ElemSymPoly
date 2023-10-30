
#include <polynomial.h>
#include <elementary_symmetric_polynomial.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <common.h>
#include <stdlib.h>

#define data_t double
#define array_t data_t *
#define matrix_t data_t **

#define state_t elementary_symmetric_state_t

struct state_t {

    int N;
    array_t elements;

    matrix_t aux_polys;

    poly_mul_state_t poly_mul_state;

};


state_t  init_elementary_symmetric_state(int size){

    if (size & (size - 1) != 0){
        size = pow(2,ceil(log2(size)));
    }

    state_t  state = (state_t ) malloc(sizeof(struct state_t));

    state->N = size;


    ALLOC_ALIGNED(state->elements, 32, sizeof(data_t) * size);

    state->aux_polys = (matrix_t) malloc(sizeof(array_t)* size);

    for (int i = 0;i < size; i++){
        ALLOC_ALIGNED(state->aux_polys[i], 32, sizeof(data_t) * 2 * size);
    }

    state->poly_mul_state = init_polynomial_mul_state(size);

    return state;

}

void update_elementary_symmetric_state(state_t state, array_t input_elements, int size){

    if (size & (size - 1) != 0){

        int current_size = size;
        size = pow(2,ceil(log2(size)));
        state->N = size;

        memcpy(state->elements, input_elements, sizeof(data_t) * current_size);
        memset(state->elements + current_size, 0, sizeof(data_t) * (size - current_size));

    } else{
        state->N = size;
        
        memcpy(state->elements, input_elements, sizeof(data_t) * state->N);

    }

    for (int i = 0;i < state->N;i++){
        state->aux_polys[i][0] = 1;
        state->aux_polys[i][1] = input_elements[i];
    }
}

static void _elem_symm_poly_comp(poly_mul_state_t poly_mul_state, matrix_t poly_matrix, int stride, int n){

    if ( n == 2){

        data_t poly_matrix_0_1 = poly_matrix[0][1];
        
        poly_matrix[0][0] = 1;
        poly_matrix[0][1] = poly_matrix_0_1 + poly_matrix[stride][1]; 
        poly_matrix[0][2] = poly_matrix_0_1 * poly_matrix[stride][1];
        poly_matrix[0][3] = 0;


    } else if (n == 4){

        _elem_symm_poly_comp(poly_mul_state, poly_matrix, 1, 2);
        _elem_symm_poly_comp(poly_mul_state, poly_matrix + 2, 1 , 2);

        data_t poly_matrix_0_1 = poly_matrix[0][1];
        data_t poly_matrix_0_2 = poly_matrix[0][2];
        data_t poly_matrix_0_3 = poly_matrix[0][3];
        
        poly_matrix[0][0] = 1;
        poly_matrix[0][1] = poly_matrix_0_1 + poly_matrix[stride][1];
        poly_matrix[0][2] = poly_matrix_0_1 * poly_matrix[stride][1] + poly_matrix_0_2 + poly_matrix[stride][2];
        //poly_matrix[0][2] = fma(poly_matrix_0_1, poly_matrix[stride][1], poly_matrix_0_2 + poly_matrix[stride][2]);     
        poly_matrix[0][3] = poly_matrix_0_1 * poly_matrix[stride][2] + poly_matrix_0_2 * poly_matrix[stride][1] + poly_matrix_0_3 + poly_matrix[stride][3];
        //poly_matrix[0][3] = fma(poly_matrix_0_1,poly_matrix[stride][2], fma(poly_matrix_0_2, poly_matrix[stride][1], poly_matrix_0_3 + poly_matrix[stride][3]));
        poly_matrix[0][4] = poly_matrix_0_1 * poly_matrix[stride][3] + poly_matrix_0_2 * poly_matrix[stride][2] + poly_matrix_0_3 * poly_matrix[stride][1];
        //poly_matrix[0][4] = fma(poly_matrix_0_1, poly_matrix[stride][3], fma(poly_matrix_0_2, poly_matrix[stride][2], poly_matrix_0_3 * poly_matrix[stride][1]));
        poly_matrix[0][5] = poly_matrix_0_2 * poly_matrix[stride][3] + poly_matrix_0_3 * poly_matrix[stride][2];
        //poly_matrix[0][5] = fma(poly_matrix_0_2, poly_matrix[stride][3], poly_matrix_0_3 * poly_matrix[stride][2]);
        poly_matrix[0][6] = poly_matrix_0_3 * poly_matrix[stride][3];
        poly_matrix[0][7] = 0;


    } else{
        _elem_symm_poly_comp(poly_mul_state, poly_matrix, stride/2, n/2);
        _elem_symm_poly_comp(poly_mul_state, poly_matrix + n/2, stride/2 , n/2);

        for (int i =0 ;i < n; i+=2*stride){

            update_polynomial_mul_state(poly_mul_state, poly_matrix[i], poly_matrix[i + stride], n );
            array_t result = polynomial_multiply(poly_mul_state);
            memcpy(poly_matrix[i], result, 2 * n * sizeof(data_t));
            
        }

    }   
}

array_t compute_elementary_symmetric_polynomials(state_t  state){

    int N = state->N;

    _elem_symm_poly_comp(state->poly_mul_state, state->aux_polys, N/2, N);

    return state->aux_polys[0];

}

void free_elementary_symmetric_state(state_t  state){
    
    free_polynomial_mul_state(state->poly_mul_state);

    for (int i = 0;i < state->N; i++){
        free(state->aux_polys[i]);
    }

    free(state->aux_polys);

    free(state->elements);

    free(state);
}


#undef state_t
#undef data_t
#undef array_t
#undef matrix_t