
#include <dft.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <panopticon.h>
#include <vecmath.h>
#include <arena.h>



typedef struct polynomial_state_d{

    int N;

    double * A_reals;
    double * A_imags;

    double * B_reals;
    double * B_imags;

    double * A_out_reals;
    double * A_out_imags;

    double * B_out_reals;
    double * B_out_imags;

    double * C_out_reals;
    double * C_out_imags;

    double ** w_reals;
    double ** w_imags;
    double ** w_reals_inverse;
    double ** w_imags_inverse;

    double * temp_C_reals;
    double * temp_C_imags;
} poly_mul_state;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"

#define ALLOC_ALIGNED(ptr, alignment, size) posix_memalign((void **)&(ptr), alignment, (size))

#pragma GCC diagnostic pop

poly_mul_state * init_polynomial_mul_state_d(int poly_size){

    int N = 2* poly_size;
    int alignment = 32;

    poly_mul_state * state = malloc(sizeof(poly_mul_state));
    
    ALLOC_ALIGNED(state->A_reals, alignment, sizeof(double) * N);
    //ALLOC_ALIGNED(state->A_imags, alignment, sizeof(double) * N);

    ALLOC_ALIGNED(state->B_reals, alignment, sizeof(double) * N);
    //ALLOC_ALIGNED(state->B_imags, alignment, sizeof(double) * N);

    ALLOC_ALIGNED(state->temp_C_reals, alignment, sizeof(double) * N);
    ALLOC_ALIGNED(state->temp_C_imags, alignment, sizeof(double) * N);
    
    ALLOC_ALIGNED(state->A_out_reals, alignment, sizeof(double) * N);
    ALLOC_ALIGNED(state->A_out_imags, alignment, sizeof(double) * N);

    ALLOC_ALIGNED(state->B_out_reals, alignment, sizeof(double) * N);
    ALLOC_ALIGNED(state->B_out_imags, alignment, sizeof(double) * N);

    ALLOC_ALIGNED(state->C_out_reals, alignment, sizeof(double) * N);
    ALLOC_ALIGNED(state->C_out_imags, alignment, sizeof(double) * N);

    memset(state->C_out_reals, 0, sizeof(double) * N);
    memset(state->C_out_imags, 0,  sizeof(double) * N);

    state->w_reals = (double **) malloc(sizeof(double *) * log2(N));
    state->w_imags = (double **) malloc(sizeof(double *) * log2(N));

    state->w_reals_inverse = (double **) malloc(sizeof(double *) * log2(N));
    state->w_imags_inverse = (double **) malloc(sizeof(double *) * log2(N));

    init_look_up_table_d(N,state->w_reals,state->w_imags);
    init_look_up_inverse_d(N,state->w_reals_inverse,state->w_imags_inverse);
    state->N = N;

    return state;
}

void update_polynomial_mul_state(poly_mul_state * state, double * A, double * B, int poly_size){

    int N = 2* poly_size;

    state->N = N;

    size_t byte_size = poly_size * sizeof(double);

    memcpy(state->A_reals, A, byte_size);
    memcpy(state->B_reals, B, byte_size);

    memset(state->A_out_reals, 0, 2 * byte_size);
    memset(state->B_out_reals, 0, 2 * byte_size);

    memset(state->A_reals + poly_size, 0, byte_size);
    memset(state->B_reals + poly_size, 0, byte_size);

    //memset(state->C_out_reals, 0, 2* byte_size);
    //memset(state->C_out_imags, 0, 2* byte_size);



}

void free_polynomial_state_d(poly_mul_state * state){

    free(state->A_reals);
    free(state->A_imags);
    free(state->B_reals);
    free(state->B_imags);
    free(state->A_out_reals);
    free(state->A_out_imags);
    free(state->B_out_reals);
    free(state->B_out_imags);
    free(state->C_out_reals);
    free(state->C_out_imags);
    free(state->w_reals);
    free(state->w_imags);
    free(state->w_reals_inverse);
    free(state->w_imags_inverse);
    free(state);

}


double * polynomial_multiply_d(poly_mul_state * state){

    int N = state->N;

    recursive_rfft_half_zero_d(state->A_reals,state->B_reals,state->A_out_reals,state->A_out_imags,state->w_reals,state->w_imags,1,N);
    //recursive_rfft_half_zero_d(state->B_reals,state->B_imags,state->B_out_reals,state->B_out_imags,state->w_reals,state->w_imags,1,N);

    for (int i = 0;i < N;i+=4){
        
        complex_4 a = LOAD_4(&state->A_out_reals[i],&state->A_out_imags[i]);
        complex_4 b = LOAD_4(&state->B_out_reals[i],&state->B_out_imags[i]);

        STORE_4(&state->temp_C_reals[i],&state->temp_C_imags[i],MUL_4(a,b));
    }

    recursive_inverse_fft_d(state->temp_C_reals,state->temp_C_imags,state->C_out_reals,state->C_out_imags,state->w_reals_inverse,state->w_imags_inverse,1,N);

    for (int i = 0;i < N;i++){
        state->C_out_reals[i] /= N;
    }

    return state->C_out_reals;

}



void multiplyPolynomials(double* poly1, int size1, double* poly2, int size2, double* result) {
    int resultSize = size1 + size2 - 1;  // Size of the resulting polynomial

    // Initialize the result polynomial with zeros
    for (int i = 0; i < resultSize; i++) {
        result[i] = 0.0;
    }

    // Perform polynomial multiplication
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            result[i + j] += poly1[i] * poly2[j];
        }
    }
}

void elementary_symmetric_recursive(poly_mul_state * state, double ** poly_array,int stride, int n){

    if ( n == 2){

        double poly_array_0_1 = poly_array[0][1];
        
        poly_array[0][0] = 1;
        poly_array[0][1] = poly_array_0_1 + poly_array[stride][1]; 
        poly_array[0][2] = poly_array_0_1 * poly_array[stride][1];
        poly_array[0][3] = 0;

    } else if (n == 4){

        elementary_symmetric_recursive(state, poly_array, 1, 2);
        elementary_symmetric_recursive(state, poly_array + 2, 1 , 2);

        double poly_array_0_1 = poly_array[0][1];
        double poly_array_0_2 = poly_array[0][2];
        double poly_array_0_3 = poly_array[0][3];

        poly_array[0][0] = 1;
        poly_array[0][1] = poly_array_0_1 + poly_array[stride][1];
        poly_array[0][2] = poly_array_0_1 * poly_array[stride][1] + poly_array_0_2 + poly_array[stride][2];
        poly_array[0][3] = poly_array_0_1 * poly_array[stride][2] + poly_array_0_2 * poly_array[stride][1] + poly_array_0_3 + poly_array[stride][3];
        poly_array[0][4] = poly_array_0_1 * poly_array[stride][3] + poly_array_0_2 * poly_array[stride][2] + poly_array_0_3 * poly_array[stride][1];
        poly_array[0][5] = poly_array_0_2 * poly_array[stride][3] + poly_array_0_3 * poly_array[stride][2];
        poly_array[0][6] = poly_array_0_3 * poly_array[stride][3];
        poly_array[0][7] = 0;

    } else{
        elementary_symmetric_recursive(state, poly_array, stride/2, n/2);
        elementary_symmetric_recursive(state, poly_array + n/2, stride/2 , n/2);

        for (int i =0 ;i < n; i+=2*stride){

            update_polynomial_mul_state(state, poly_array[i], poly_array[i + stride], n / 2);
            double * result = polynomial_multiply_d(state);
            memcpy(poly_array[i], result, n / 2 * sizeof(double));
            
        }

    }   

}


int main(){
    GLOBAL_TIMER(MICROSECONDS,CLOCK_MONOTONIC_RAW)

    

    int N = pow(2,14);
/*
    p->degree_bound = 2*N;
    a->degree_bound = 2*N;

    posix_memalign((void**)&p->coefficients, 32, 2 *N * sizeof(double));
    posix_memalign((void **)&a->coefficients, 32, 2 *N * sizeof(float));

    for (int i = 0;i < N;i++){
        p->coefficients[i] = i;
        a->coefficients[i] = i;
    }

    struct polynomial_d * q = malloc(sizeof(struct polynomial_d));
    struct polynomial * b = malloc(sizeof(struct polynomial));

    q->degree_bound = 2*N;
    b->degree_bound = 2*N;
    posix_memalign((void**)&q->coefficients, 32, 2*N * sizeof(double));
    posix_memalign((void **)&b->coefficients, 32, 2*N *sizeof(float));


    for (int i = 0;i < N;i++){
        q->coefficients[i] = i;
        b->coefficients[i] = i;
    }

    struct polynomial_d * r = malloc(sizeof(struct polynomial_d));
    struct polynomial * c = malloc(sizeof(struct polynomial));

    r->degree_bound = 2*N;
    c->degree_bound = 2*N;
    posix_memalign((void**)&r->coefficients, 32, 2*N * sizeof(double));
    posix_memalign((void **)&c->coefficients, 32, 2*N * sizeof(float));


    WATCH("fft-poly-init")
    struct polynomial_state_d * state = init_polynomial_mul_state_d(p,q,r);
    STOP_WATCH("fft-poly-init")
    WATCH("fft-poly-init-f")
    struct polynomial_state * state_f = init_polynomial_mul_state(a,b,c);
    STOP_WATCH("fft-poly-init-f")
    WATCH("fft-poly")
    polynomial_multiply_d(state);

    STOP_WATCH("fft-poly")

    WATCH("fft-poly-f")
    polynomial_multiply(state_f);
    STOP_WATCH("fft-poly-f")

    WATCH("poly-naive")

    //multiplyPolynomials(p->coefficients,N,q->coefficients,N,r->coefficients);

    STOP_WATCH("poly-naive")
*/
    int aligement = 32;
    double * elements;
    double * elements_1;
    volatile double * out;

    ALLOC_ALIGNED(elements, aligement, 2*N*sizeof(double));
    ALLOC_ALIGNED(elements_1, aligement, 2*N*sizeof(double));
    ALLOC_ALIGNED(out, aligement, 2*N* sizeof(double));

    double ** polys = (double **)malloc(sizeof(double *)*N);

    for (int i = 0;i < N;i++){
        polys[i];
        ALLOC_ALIGNED(polys[i],aligement,sizeof(double)* 2* N);
        elements[i] = i;
        elements_1[i] = i;
        polys[i][0] = 1;
        polys[i][1] = 1;
        out[i] =0;
    }

    WATCH("poly-elem")
    struct polynomial_state_d * state = init_polynomial_mul_state_d(N);
    
    //elementary_symmetric_recursive(arena,polys,N/2,N);
    STOP_WATCH("poly-elem")
/*
    for (int i =0;i< 2*N;i++){
        printf("Val: %lf \n",polys[0][i]);
    }
*/  
    WATCH("poly_mul")
    for (int i =0;i< 10;i++){
        update_polynomial_mul_state(state,elements,elements_1, N);
        double * ret =polynomial_multiply_d(state);
    }
    
    STOP_WATCH("poly_mul");



//    for (int i =0;i< 2;i++){
//        printf("Val: %lf  \n",ret[i]);
//    }


    
    WATCH("naive")
    for (int i =0;i< 10;i++){
        multiplyPolynomials(elements,N,elements_1,N,out);
    }
    STOP_WATCH("naive")

    return 0;

}

