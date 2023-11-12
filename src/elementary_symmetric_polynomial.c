
#include <polynomial.h>
#include <elementary_symmetric_polynomial.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <common.h>
#include <stdlib.h>


# define data_t mpfr_t
#define array_t data_t *
#define matrix_t data_t **

#define state_t elementary_symmetric_state_t

struct state_t {

    int N;
    double * elements;

    matrix_t aux_polys;

    poly_mul_state_t poly_mul_state;

};


state_t  init_elementary_symmetric_state(int size){

    if (size & (size - 1) != 0){
        size = pow(2,ceil(log2(size)));
    }

    state_t  state =  ALLOC(sizeof(struct state_t));

    state->N = size;


    state->elements = ALLOC(sizeof(double) * size);

    state->aux_polys = ALLOC(sizeof(array_t) * size);

    for (int i = 0;i < size; i++){
        state->aux_polys[i] = ALLOC(sizeof(data_t) * 2 * size);
        for (int j =0; j < size *2; j++){
            mpfr_init2(state->aux_polys[i][j],PRECISION);
        }
    }

    state->poly_mul_state = init_polynomial_mul_state(size);

    return state;

}

void update_elementary_symmetric_state(state_t state, double * input_elements, int size){

    if (size & (size - 1) != 0){

        int current_size = size;
        size = pow(2,ceil(log2(size)));
        state->N = size;

        memcpy(state->elements, input_elements, sizeof(double) * current_size);
        memset(state->elements + current_size, 0, sizeof(double) * (size - current_size));

    } else{
        state->N = size;
        
        memcpy(state->elements, input_elements, sizeof(double ) * state->N);

    }

    for (int i = 0;i < state->N;i++){
		//mpfr_init2(state->aux_polys[i][0],PRECISION);
		//mpfr_init2(state->aux_polys[i][1],PRECISION);
		mpfr_set_d(state->aux_polys[i][0],1, MPFR_RNDN);
		mpfr_set_d(state->aux_polys[i][1],input_elements[i], MPFR_RNDN);
		//mpfr_printf("%Rf\n", state->aux_polys[i][1]);
    }
}

static void elem_symm_poly_comp(poly_mul_state_t poly_mul_state, matrix_t poly_matrix, int stride, int n){

    if ( n == 2){


        //data_t poly_matrix_0_1 = poly_matrix[0][1];
        //poly_matrix[0][1] = poly_matrix_0_1 + poly_matrix[stride][1];
        //poly_matrix[0][2] = poly_matrix_0_1 * poly_matrix[stride][1];
		mpfr_set_d(poly_matrix[0][0],1, MPFR_RNDN);
		mpfr_mul(poly_matrix[0][2],poly_matrix[0][1],poly_matrix[stride][1], MPFR_RNDN);
		mpfr_add(poly_matrix[0][1],poly_matrix[0][1],poly_matrix[stride][1],MPFR_RNDN);
		mpfr_set_d(poly_matrix[0][3],0,MPFR_RNDN);

	} else if (n == 4){

        elem_symm_poly_comp(poly_mul_state, poly_matrix, 1, 2);
        elem_symm_poly_comp(poly_mul_state, poly_matrix + 2, 1 , 2);

        //data_t poly_matrix_0_1 = poly_matrix[0][1];
        //data_t poly_matrix_0_2 = poly_matrix[0][2];
        //data_t poly_matrix_0_3 = poly_matrix[0][3];

		data_t poly_matrix_0_1, poly_matrix_0_2, poly_matrix_0_3, aux_0, aux_1, aux_2;
		mpfr_inits2(PRECISION, poly_matrix_0_1, poly_matrix_0_2, poly_matrix_0_3,aux_0,aux_1,aux_2, (mpfr_ptr)NULL);
		mpfr_set(poly_matrix_0_1, poly_matrix[0][1], MPFR_RNDN);
		mpfr_set(poly_matrix_0_2, poly_matrix[0][2], MPFR_RNDN);
		mpfr_set(poly_matrix_0_3, poly_matrix[0][3], MPFR_RNDN);

		mpfr_set_d(poly_matrix[0][0],1, MPFR_RNDN);

		mpfr_add(poly_matrix[0][1],poly_matrix_0_1,poly_matrix[stride][1],MPFR_RNDN);

		mpfr_mul(aux_0,poly_matrix_0_1,poly_matrix[stride][1],MPFR_RNDN);
		mpfr_add(aux_1, poly_matrix_0_2, poly_matrix[stride][2], MPFR_RNDN);
		mpfr_add(poly_matrix[0][2], aux_0, aux_1, MPFR_RNDN);

		mpfr_mul(aux_0,poly_matrix_0_1,poly_matrix[stride][2],MPFR_RNDN);
		mpfr_mul(aux_1,poly_matrix_0_2,poly_matrix[stride][1],MPFR_RNDN);
		mpfr_add(aux_2, poly_matrix_0_3, poly_matrix[stride][3], MPFR_RNDN);
		mpfr_add(poly_matrix[0][3], aux_0, aux_1, MPFR_RNDN);
		mpfr_add(poly_matrix[0][3], poly_matrix[0][3], aux_2, MPFR_RNDN);

		mpfr_mul(aux_0,poly_matrix_0_1,poly_matrix[stride][3],MPFR_RNDN);
		mpfr_mul(aux_1,poly_matrix_0_2,poly_matrix[stride][2],MPFR_RNDN);
		mpfr_mul(aux_2,poly_matrix_0_3,poly_matrix[stride][1],MPFR_RNDN);
		mpfr_add(poly_matrix[0][4], aux_0, aux_1, MPFR_RNDN);
		mpfr_add(poly_matrix[0][4], poly_matrix[0][4], aux_2, MPFR_RNDN);

		mpfr_mul(aux_0,poly_matrix_0_2,poly_matrix[stride][3],MPFR_RNDN);
		mpfr_mul(aux_1,poly_matrix_0_3,poly_matrix[stride][2],MPFR_RNDN);
		mpfr_add(poly_matrix[0][5], aux_0, aux_1, MPFR_RNDN);

		mpfr_mul(poly_matrix[0][6],poly_matrix_0_3,poly_matrix[stride][3],MPFR_RNDN);
		mpfr_set_d(poly_matrix[0][7],0,MPFR_RNDN);

		mpfr_clears(poly_matrix_0_1, poly_matrix_0_2, poly_matrix_0_3,aux_0,aux_1,aux_2, (mpfr_ptr)NULL);
        //poly_matrix[0][1] = poly_matrix_0_1 + poly_matrix[stride][1];
        //poly_matrix[0][2] = poly_matrix_0_1 * poly_matrix[stride][1] + poly_matrix_0_2 + poly_matrix[stride][2];
        //poly_matrix[0][3] = poly_matrix_0_1 * poly_matrix[stride][2] + poly_matrix_0_2 * poly_matrix[stride][1] + poly_matrix_0_3 + poly_matrix[stride][3];
        //poly_matrix[0][4] = poly_matrix_0_1 * poly_matrix[stride][3] + poly_matrix_0_2 * poly_matrix[stride][2] + poly_matrix_0_3 * poly_matrix[stride][1];
        //poly_matrix[0][5] = poly_matrix_0_2 * poly_matrix[stride][3] + poly_matrix_0_3 * poly_matrix[stride][2];
        //poly_matrix[0][6] = poly_matrix_0_3 * poly_matrix[stride][3];


    } else{
        elem_symm_poly_comp(poly_mul_state, poly_matrix, stride/2, n/2);
        elem_symm_poly_comp(poly_mul_state, poly_matrix + n/2, stride/2 , n/2);

        for (int i =0 ;i < n; i+=2*stride){

            update_polynomial_mul_state(poly_mul_state, poly_matrix[i], poly_matrix[i + stride], n );
            array_t result = polynomial_multiply(poly_mul_state);
			for (int j = 0;j < 2*n;j++){
				mpfr_set(poly_matrix[i][j], result[j], MPFR_RNDN);
			}

            
        }

    }   
}


double * compute_elementary_symmetric_polynomials(state_t  state){

    int N = state->N;

	elem_symm_poly_comp(state->poly_mul_state, state->aux_polys, N/2, N);

	double * result = ALLOC(sizeof(double) * (N+1));
	for (int i = 0;i < N+1;i++){
		//mpfr_printf("%Rf\n", state->aux_polys[0][i]);
		result[i] = mpfr_get_d(state->aux_polys[0][i], MPFR_RNDN);
	}
    return result;

}

void free_elementary_symmetric_state(state_t  state){
    
    free_polynomial_mul_state(state->poly_mul_state);

    int size = state->N;

    for (int i = 0;i < size; i++){
        for (int j =0; j < 2*size;j++){
            mpfr_clear(state->aux_polys[i][j]);
        }
        free(state->aux_polys[i]);
    }

    free(state->aux_polys);

    free(state->elements);

    free(state);

    mpfr_free_cache();
}


#undef state_t
#undef data_t
#undef array_t
#undef matrix_t