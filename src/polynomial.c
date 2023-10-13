
#include <dft.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <panopticon.h>
#include <vecmath.h>
#include <arena.h>

struct polynomial_d{
    int degree_bound;
    double * coefficients;
};

struct polynomial{
    int degree_bound;
    float * coefficients;
};

struct polynomial_state_d{

    int N;

    double * A_reals;
    double * A_imags;

    double * B_reals;
    double * B_imags;

    double * C_reals;
    double * C_imags;

    double * A_out_reals;
    double * A_out_imags;

    double * B_out_reals;
    double * B_out_imags;

    double * C_out_reals;
    double * C_out_imags;

    double * w_reals;
    double * w_imags;
    double * w_reals_inverse;
    double * w_imags_inverse;
};

struct polynomial_state{

    int N;

    float * A_reals;
    float * A_imags;

    float * B_reals;
    float * B_imags;

    float * C_reals;
    float * C_imags;

    float * A_out_reals;
    float * A_out_imags;

    float * B_out_reals;
    float * B_out_imags;

    float * C_out_reals;
    float * C_out_imags;

    float * w_reals;
    float * w_imags;
    float * w_reals_inverse;
    float * w_imags_inverse;

};

struct polynomial_state_d * init_polynomial_mul_state_d(struct polynomial_d * A, struct polynomial_d * B, struct polynomial_d * C){

    int N = A->degree_bound;

    struct polynomial_state_d * state = malloc(sizeof(struct polynomial_state_d));

    state->N = N;
    posix_memalign((void**)&state->A_reals, 32, N * sizeof(double));
    posix_memalign((void**)&state->A_imags, 32, N * sizeof(double));
    posix_memalign((void**)&state->B_reals, 32, N * sizeof(double));
    posix_memalign((void**)&state->B_imags, 32, N * sizeof(double));
    posix_memalign((void**)&state->C_reals, 32, N * sizeof(double));
    posix_memalign((void**)&state->C_imags, 32, N * sizeof(double));
    posix_memalign((void**)&state->A_out_reals, 32, N * sizeof(double));
    posix_memalign((void**)&state->A_out_imags, 32, N * sizeof(double));
    posix_memalign((void**)&state->B_out_reals, 32, N * sizeof(double));
    posix_memalign((void**)&state->B_out_imags, 32, N * sizeof(double));
    posix_memalign((void**)&state->C_out_reals, 32, N * sizeof(double));
    posix_memalign((void**)&state->C_out_imags, 32, N * sizeof(double));
    posix_memalign((void**)&state->w_reals, 32, N * sizeof(double));
    posix_memalign((void**)&state->w_imags, 32, N * sizeof(double));
    posix_memalign((void**)&state->w_reals_inverse, 32, N * sizeof(double));
    posix_memalign((void**)&state->w_imags_inverse, 32, N * sizeof(double));

    memset(state->A_imags,0,N*sizeof(double));
    memset(state->B_imags,0,N*sizeof(double));

    memcpy(state->A_reals,A->coefficients,N*sizeof(double));
    memcpy(state->B_reals,B->coefficients,N*sizeof(double));


    init_look_up_table_d(N,state->w_reals,state->w_imags);
    init_look_up_inverse_d(N,state->w_reals_inverse,state->w_imags_inverse);

    return state;
}

void update_polynomial_mul_state_d(struct polynomial_state_d * state, struct polynomial_d * A, struct polynomial_d * B, struct polynomial_d * C){


}


struct polynomial_state * init_polynomial_mul_state(struct polynomial * A, struct polynomial * B, struct polynomial * C){

    int N = A->degree_bound;

    struct polynomial_state * state = malloc(sizeof(struct polynomial_state));

    state->N = N;
    posix_memalign((void**)&state->A_reals, 32, N * sizeof(float));
    posix_memalign((void**)&state->A_imags, 32, N * sizeof(float));
    posix_memalign((void**)&state->B_reals, 32, N * sizeof(float));
    posix_memalign((void**)&state->B_imags, 32, N * sizeof(float));
    posix_memalign((void**)&state->C_reals, 32, N * sizeof(float));
    posix_memalign((void**)&state->C_imags, 32, N * sizeof(float));
    posix_memalign((void**)&state->A_out_reals, 32, N * sizeof(float));
    posix_memalign((void**)&state->A_out_imags, 32, N * sizeof(float));
    posix_memalign((void**)&state->B_out_reals, 32, N * sizeof(float));
    posix_memalign((void**)&state->B_out_imags, 32, N * sizeof(float));
    posix_memalign((void**)&state->C_out_reals, 32, N * sizeof(float));
    posix_memalign((void**)&state->C_out_imags, 32, N * sizeof(float));
    posix_memalign((void**)&state->w_reals, 32, N * sizeof(float));
    posix_memalign((void**)&state->w_imags, 32, N * sizeof(float));
    posix_memalign((void**)&state->w_reals_inverse, 32, N * sizeof(float));
    posix_memalign((void**)&state->w_imags_inverse, 32, N * sizeof(float));

    memset(state->A_imags,0,N*sizeof(float));
    memset(state->B_imags,0,N*sizeof(float));

    memset(state->A_out_reals,0,N*sizeof(float));
    memset(state->A_out_imags,0,N*sizeof(float));

    memcpy(state->A_reals,A->coefficients,N*sizeof(float));
    memcpy(state->B_reals,B->coefficients,N*sizeof(float));


    init_look_up_table(N,state->w_reals,state->w_imags);
    init_look_up_inverse(N,state->w_reals_inverse,state->w_imags_inverse);

    return state;

}

void free_polynomial_state_d(struct polynomial_state_d * state){

    free(state->A_reals);
    free(state->A_imags);
    free(state->B_reals);
    free(state->B_imags);
    free(state->C_reals);
    free(state->C_imags);
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


void polynomial_multiply_d(struct polynomial_state_d * state){

    int N = state->N;

    recursive_fft_d(state->A_reals,state->A_imags,state->A_out_reals,state->A_out_imags,state->w_reals,state->w_imags,1,N);
    recursive_fft_d(state->B_reals,state->B_imags,state->B_out_reals,state->B_out_imags,state->w_reals,state->w_imags,1,N);

    for (int i = 0;i < N;i+=4){
        
        complex_4 a = LOAD_4(&state->A_out_reals[i],&state->A_out_imags[i]);
        complex_4 b = LOAD_4(&state->B_out_reals[i],&state->B_out_imags[i]);

        STORE_4(&state->C_reals[i],&state->C_imags[i],MUL_4(a,b));
    }

    recursive_inverse_fft_d(state->C_reals,state->C_imags,state->C_out_reals,state->C_out_imags,state->w_reals_inverse,state->w_imags_inverse,1,N);

    //state->C_reals = state->C_out_reals;
    memcpy(state->C_reals,state->C_out_reals,N*sizeof(double));
    for (int i = 0;i < N;i++){
        state->C_reals[i] /= N;
    }


    

}

void polynomial_multiply(struct polynomial_state * state){

    int N = state->N;



    recursive_fft(state->A_reals,state->A_imags,state->A_out_reals,state->A_out_imags,state->w_reals,state->w_imags,1,N);
    recursive_fft(state->B_reals,state->B_imags,state->B_out_reals,state->B_out_imags,state->w_reals,state->w_imags,1,N);

    for (int i = 0;i < N;i+=8){
        
        complex_8 a = LOAD(&state->A_out_reals[i],&state->A_out_imags[i]);
        complex_8 b = LOAD(&state->B_out_reals[i],&state->B_out_imags[i]);

        STORE(&state->C_reals[i],&state->C_imags[i],MUL(a,b));
    }

    recursive_inverse_fft(state->C_reals,state->C_imags,state->C_out_reals,state->C_out_imags,state->w_reals_inverse,state->w_imags_inverse,1,N);

    for (int i = 0;i < N;i++){
        state->C_out_reals[i] /= N;
    }

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

void elementary_symmetric_recursive(struct polynomial_d ** poly_array, int stride, int n){

    if (stride > n/2){
        return;
    } else {
        elementary_symmetric_recursive(poly_array,stride << 1, n);

        for (int i = 0; i < n; i+=2*stride ){
            struct polynomial_state_d * state = init_polynomial_mul_state_d(poly_array[i],poly_array[i+stride],poly_array[i]);
            polynomial_multiply_d(state);
        }
    }


}


int main(){
    GLOBAL_TIMER(MICROSECONDS,CLOCK_MONOTONIC_RAW)


    struct polynomial_d * p = malloc(sizeof(struct polynomial_d));
    struct polynomial * a = malloc(sizeof(struct polynomial));

    int N = pow(2,6);

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

    double * elements;
    posix_memalign((void**)&elements, 32, N * sizeof(double));

    struct polynomial_d ** polys = malloc(sizeof(struct polynomial_d*)*N);

    for (int i = 0;i < N;i++){
        struct polynomial_d * poly = malloc(sizeof(struct polynomial_d));
        posix_memalign((void**)&poly->coefficients, 32, N * sizeof(double));
        poly->degree_bound = N;
        elements[i] = i;
        poly->coefficients[0] = 1;
        poly->coefficients[1] = i; 
        polys[i] = poly;
    }


    WATCH("poly-elem")
    elementary_symmetric_recursive(polys,1,N);
    STOP_WATCH("poly-elem")


    free_polynomial_state_d(state);
    free(p->coefficients);
    free(p);
    free(q->coefficients);
    free(q);
    free(r->coefficients);
    free(r);

    return 0;

}

