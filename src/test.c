#include <stdio.h>
#include <time.h>
#include <cblas.h>
#include <lapacke.h>
#include <dichotomus.h>
#include <utils.h>
#include <panopticon.h>
#include <roots.h>
#include <dft.h>
#include <vecmath.h>
#include <arena.h>

#define W(N,k) (cexp(-2.0f * M_PI * I * (float)(k)/(float)(N)))


typedef struct complex_array{

    float * reals;
    float * imags;

} complex_array;

static inline complex_array * shift_array(complex_array * array, int n){
    array->imags = array->imags + n;
    array->reals = array->reals + n;
    return array;
}

complex_array * get_even(complex_array * A,Arena_T arena, int n){

    //complex_array * ret = malloc(sizeof(complex_array));
    complex_array * ret = arena_alloc(arena,sizeof(complex_array));
    //float * reals = calloc(n/2,sizeof(float));
    //float * imags = calloc(n/2,sizeof(float));
    float * reals = arena_alloc(arena,sizeof(float)*n/2);
    float * imags = arena_alloc(arena,sizeof(float)*n/2);

    ret->reals = reals;
    ret->imags = imags;

    for (int i =0;i < n/2;i++){
        ret->reals[i] = A->reals[2*i];
        ret->imags[i] = A->imags[2*i];
    }
    

    

    return ret;
}

complex_array * get_odd(complex_array * A, Arena_T arena, int n){

    //complex_array * ret = malloc(sizeof(complex_array));
    //float * reals = calloc(n/2,sizeof(float));
    //float * imags = calloc(n/2,sizeof(float));
    complex_array * ret = arena_alloc(arena,sizeof(complex_array));
    float * reals = arena_alloc(arena,sizeof(float)*n/2);
    float * imags = arena_alloc(arena,sizeof(float)*n/2);
    ret->reals = reals;
    ret->imags = imags;

    for (int i =0;i < n/2;i++){
        ret->reals[i] = A->reals[2*i+1];
        ret->imags[i] = A->imags[2*i+1];
    }

    return ret;
}



void recursive_fft( float * restrict in_reals,  float *  restrict in_imags,  float * restrict out_reals,  float * restrict out_imags, complex_array * restrict look_up_table, int stride, int n){


    if (n ==2){


        float o1_real = in_reals[0];
        float o1_imag = in_imags[0];
        float o2_real = in_reals[stride];
        float o2_imag = in_imags[stride];

        out_reals[0] = o1_real + o2_real;
        out_imags[0] = o1_imag + o2_imag;
        out_reals[n/2] = o1_real - o2_real;
        out_imags[n/2] = o1_imag - o2_imag;
        
    } else if ( n == 4){

        recursive_fft(in_reals,in_imags, out_reals, out_imags, look_up_table, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, look_up_table, stride << 1, n >> 1);
        

        float o0_real = out_reals[0];
        float o0_imag = out_imags[0];
        float o1_real = out_reals[1];
        float o1_imag = out_imags[1];
        float o2_real = out_reals[2];
        float o2_imag = out_imags[2];
        float o3_real = out_reals[3];
        float o3_imag = out_imags[3];
        
        

        out_reals[0] = o0_real + o2_real;
        out_imags[0] = o0_imag + o2_imag;

        out_reals[1] = o1_real - o3_imag;
        out_imags[1] = o1_imag + o3_real;

        out_reals[2] = o0_real - o2_real;
        out_imags[2] = o0_imag - o2_imag;

        out_reals[3] = o1_real + o3_imag;
        out_imags[3] = o1_imag - o3_real;



    } else if ( n == 8){
        
        float w_n_real = cos(2.0f * M_PI * (float)(1)/(float)(8));
        float w_n_imag = sin(2.0f * M_PI * (float)(1)/(float)(8));
        float w_real =1;
        float w_imag =0;

        recursive_fft(in_reals,in_imags, out_reals, out_imags, look_up_table, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, look_up_table, stride << 1, n >> 1);

        for (int k =0; k < n/2; k++){
            //printf("Val: %f + i%f \n",out->reals[k],out->imags[k]);
            float t_real = w_real*out_reals[k + n/2] - w_imag*out_imags[k + n/2];
            float t_imag = w_real*out_imags[k + n/2] + w_imag*out_reals[k + n/2];

            float temp_real = out_reals[k];
            float temp_imag = out_imags[k];
            
            out_reals[k] = temp_real + t_real;
            out_imags[k] = temp_imag + t_imag;
            out_reals[k+n/2] = temp_real - t_real;
            out_imags[k+n/2] = temp_imag - t_imag;


            float prev_w_real = w_real;
            float prev_w_imag = w_imag;

            // Update w_real and w_imag using the correct formulas
            w_real = prev_w_real * w_n_real - prev_w_imag * w_n_imag;
            w_imag = prev_w_real * w_n_imag + prev_w_imag * w_n_real;
        }
        
        
    } else if ( n == 16){
        recursive_fft(in_reals,in_imags, out_reals, out_imags, look_up_table, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, look_up_table, stride << 1, n >> 1);

        complex_8 w, y_1_k, t, y_0_k;
        int k;

        w = LOAD(&look_up_table->reals[0],&look_up_table->imags[0]);
        y_1_k = LOAD(&out_reals[8],&out_imags[8]);
        y_0_k = LOAD(&out_reals[0],&out_imags[0]);
        
        t = MUL(w,y_1_k);
        
        STORE(&out_reals[0],&out_imags[0],ADD(y_0_k,t));
        STORE(&out_reals[8],&out_imags[8],SUB(y_0_k,t));

    } else {

        recursive_fft(in_reals,in_imags, out_reals, out_imags, look_up_table, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, look_up_table, stride << 1, n >> 1);

        complex_8 w, y_1_k, t, y_0_k;
        int k;
        
        
        for (k =0; k < n/2; k+=8){
            
            w = LOAD(&look_up_table->reals[k],&look_up_table->imags[k]);
            y_1_k = LOAD(&out_reals[k+n/2],&out_imags[k + n/2]);
            y_0_k = LOAD(&out_reals[k],&out_imags[k]);
            
            t = MUL(w,y_1_k);
            
            STORE(&out_reals[k],&out_imags[k],ADD(y_0_k,t));
            STORE(&out_reals[k+n/2],&out_imags[k+n/2],SUB(y_0_k,t));
            
        }
       
        
        
    }

}



complex_array * init_look_up_table(int N, Arena_T arena){

    //complex_array * look_up_table = malloc(sizeof(complex_array));
    //look_up_table->reals = malloc(sizeof(float)*N);
    //look_up_table->imags = malloc(sizeof(float)*N);

    complex_array * look_up_table = malloc(sizeof(complex_array));
    float * reals;
    float * imags;

    int result = posix_memalign((void**)&reals, 32, N * sizeof(float));
    int result2 = posix_memalign((void**)&imags, 32, N * sizeof(float));
    
    look_up_table->reals = reals;
    look_up_table->imags = imags;

    for (int i =0;i < N/2;i++){
        look_up_table->reals[i] = cos(2.0f * M_PI * (float)(i)/(float)(N));
        look_up_table->imags[i] = sin(2.0f * M_PI * (float)(i)/(float)(N));
    }

    return look_up_table;
}

int main(){

    GLOBAL_TIMER(MICROSECONDS,CLOCK_MONOTONIC_RAW)  
    Arena_T arena = arena_new();
    int N =pow(2,16);
    int alignment = 32;
    //complex_array * in = malloc(sizeof( complex_array ));

    //in->reals = malloc(sizeof(float)*N);
    //in->imags = calloc(N,sizeof(float));
    //complex_array * in = arena_alloc(arena,sizeof(complex_array));

    //float * in_reals = arena_calloc(arena,N,sizeof(float));
    //float *  in_imags = arena_calloc(arena,N,sizeof(float));
    float * in_reals;
    float * in_imags;
    int result = posix_memalign((void**)&in_reals, alignment, N * sizeof(float));
    int result2 = posix_memalign((void**)&in_imags, alignment, N * sizeof(float));

    //complex_array * out = malloc(sizeof( complex_array ));
    //float *  out_reals = arena_calloc(arena,N,sizeof(float));
    //float *  out_imags = arena_calloc(arena,N,sizeof(float));

    float * out_reals;
    float * out_imags;

    int result3 = posix_memalign((void**)&out_reals, alignment, N * sizeof(float));
    int result4 = posix_memalign((void**)&out_imags, alignment, N * sizeof(float));


    for (int i =0;i < N;i++){
        
        in_reals[i]= rand()% 100 +1;
    }
  
/*
    in_reals[0] = 1;
    in_reals[1] =2;
    in_reals[2] = 3;
    in_reals[3] = -4;
    in_reals[4] = 1;
    in_reals[5] =2;
    in_reals[6] = 3;
    in_reals[7] = 4;
    in_reals[8] = 1;
    in_reals[9] =2;
    in_reals[10] = 3;
    in_reals[11] = -4;
    in_reals[12] = 1;
    in_reals[13] =2;
    in_reals[14] = 3;
    in_reals[15] = 4;

 */

    WATCH("copy")
    complex_array * look_up_table = init_look_up_table(N,arena);
    STOP_WATCH("copy")
    //A = reverse_copy(vec);

    
    WATCH("dft")
    //dft(vec,A);
    recursive_fft(in_reals,in_imags,out_reals,out_imags,look_up_table,1,N);

    STOP_WATCH("dft")


    for (int i =0;i<2;i++){
        printf("Valw: %f + i%f \n",out_reals[i],out_imags[i]);
    }
    //print_complex_vector(A);
    arena_dispose(&arena);
    return 0;

}