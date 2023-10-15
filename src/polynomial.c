
#include <dft.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <panopticon.h>
#include <vecmath.h>
#include <arena.h>



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

struct polynomial_state_d * init_polynomial_mul_state_d(Arena_T arena, double * A, double * B, double * C, int poly_size){

    int N = 2* poly_size;

    
    

    /* 

    struct polynomial_state_d * state = malloc(sizeof(struct polynomial_state_d));
    
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
   */
    WATCH("fft-poly-init")
    struct polynomial_state_d * state = arena_alloc(arena,sizeof(struct polynomial_state_d));
    state->A_reals = arena_alloc(arena,N * sizeof(double));
    state->A_imags = arena_alloc(arena,N * sizeof(double));

    state->B_reals = arena_alloc(arena,N * sizeof(double));
    state->B_imags = arena_alloc(arena,N * sizeof(double));

    state->C_reals = C;
    state->C_imags = arena_alloc(arena,N * sizeof(double));

    state->A_out_reals = arena_alloc(arena,N * sizeof(double));
    state->A_out_imags = arena_alloc(arena,N * sizeof(double));

    state->B_out_reals = arena_alloc(arena,N * sizeof(double));
    state->B_out_imags = arena_alloc(arena,N * sizeof(double));

    state->C_out_reals = arena_alloc(arena,N * sizeof(double));
    state->C_out_imags = arena_alloc(arena,N * sizeof(double));

    state->w_reals = arena_alloc(arena,N * sizeof(double));
    state->w_imags = arena_alloc(arena,N * sizeof(double));

    state->w_reals_inverse = arena_alloc(arena,N * sizeof(double));
    state->w_imags_inverse = arena_alloc(arena,N * sizeof(double));

    memset(state->A_imags,0,N*sizeof(double));
    memset(state->B_imags,0,N*sizeof(double));

    memcpy(state->A_reals,A,N*sizeof(double));
    memcpy(state->B_reals,B,N*sizeof(double));
    

    init_look_up_table_d(N,state->w_reals,state->w_imags);
    init_look_up_inverse_d(N,state->w_reals_inverse,state->w_imags_inverse);
    state->N = N;
    STOP_WATCH("fft-poly-init")
    return state;
}


struct polynomial_state * init_polynomial_mul_state(float * A, float * B, float * C, int poly_size){

    int N = 2* poly_size;

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

    memcpy(state->A_reals,A,N*sizeof(float));
    memcpy(state->B_reals,B,N*sizeof(float));


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

void elementary_symmetric_recursive(Arena_T arena, double ** poly_array,int stride, int n){

    if ( n == 2){
        
        printf(" n stride : %d  %d\n", n, stride);
        puts("-----------------------------------------");
        for (int i =0;i< n;i++){
            printf("Val A: %lf ------- %lf \n",poly_array[0][i],poly_array[stride][i]);
        }
        double poly_array_0_1 = poly_array[0][1];
        
        poly_array[0][0] = 1;
        poly_array[0][1] = poly_array_0_1 + poly_array[stride][1]; 
        poly_array[0][2] = poly_array_0_1 * poly_array[stride][1];
        poly_array[0][3] = 0;
        
        puts("-----------------------------------------");
        //base
    } else if (n == 4){
        printf(" n stride : %d  %d\n", n, stride);
        puts("-----------------------------------------");

        

        elementary_symmetric_recursive(arena, poly_array, 1, 2);
        elementary_symmetric_recursive(arena, poly_array + 2, 1 , 2);

        for (int i =0;i< n;i++){
            printf("Val A: %lf ------- %lf \n",poly_array[0][i],poly_array[stride][i]);
        }

        double poly_array_0_1 = poly_array[0][1];
        double poly_array_0_2 = poly_array[0][2];
        double poly_array_0_3 = poly_array[0][3];

        poly_array[0][0] = 1;
        poly_array[0][1] = poly_array_0_1 + poly_array[stride][1];
        poly_array[0][2] = poly_array_0_1 * poly_array[stride][1] + poly_array[stride][2] + poly_array_0_2;
        poly_array[0][3] = poly_array_0_1 * poly_array[stride][2] + poly_array[stride][1] * poly_array_0_2 + poly_array_0_3 + poly_array[stride][3];
        poly_array[0][4] = poly_array_0_3 * poly_array[stride][1] + poly_array[stride][2] * poly_array_0_2 + poly_array_0_1 * poly_array[stride][3];
        poly_array[0][5] = poly_array_0_3 * poly_array[stride][2] + poly_array[stride][3] * poly_array_0_2;
        poly_array[0][6] = poly_array_0_3 * poly_array[stride][3];
        poly_array[0][7] = 0;

        puts("-----------------------------------------");

    } else{
        elementary_symmetric_recursive(arena, poly_array, stride/2, n/2);
        elementary_symmetric_recursive(arena, poly_array + n/2, stride/2 , n/2);

        printf(" n stride : %d  %d\n", n, stride);
        puts("-----------------------------------------");
        for (int i =0;i< n;i++){
            printf("Val A: %lf ------- %lf \n",poly_array[0][i],poly_array[stride][i]);
        }

        for (int i =0 ;i < n; i+=2*stride){
            printf("i and stride: %d %d\n",i,stride);
            for (int j =0; j< 4;j++){
                printf("@ %lf <-> %lf \n",poly_array[i][j],poly_array[i + stride][j]);
            }

            struct polynomial_state_d * state = init_polynomial_mul_state_d(arena, poly_array[i], poly_array[i + stride], poly_array[i], n/2);
            polynomial_multiply_d(state);

        }
        puts("-----------------------------------------");
    }   

    
    
            
    //struct polynomial_state_d * state = init_polynomial_mul_state_d(arena, poly_array[i],poly_array[i+stride],poly_array[i],n);
    //polynomial_multiply_d(state);
    


}


int main(){
    GLOBAL_TIMER(MICROSECONDS,CLOCK_MONOTONIC_RAW)


    

    int N = pow(2,3);
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
    Arena_T arena = arena_new();
    double * elements = arena_alloc(arena,N*sizeof(double));

    double ** polys = arena_alloc(arena,sizeof(double *)*N);

    for (int i = 0;i < N;i++){
        polys[i] = arena_alloc(arena,N*sizeof(double));
        elements[i] = 1;
        polys[i][0] = 1;
        polys[i][1] = 1;
    }

    
    WATCH("poly-elem")
    elementary_symmetric_recursive(arena,polys,N/2,N);
    STOP_WATCH("poly-elem")

    for (int i =0;i< N;i++){
        printf("Val: %lf \n",polys[0][i]);
    }

    return 0;

}

