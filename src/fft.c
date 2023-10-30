#include <fft.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vecmath.h>
#include <common.h>
#include <panopticon.h> 

#define data_t double
#define array_t data_t * restrict
#define matrix_t data_t ** restrict


void init_look_up_table_d(int N, matrix_t reals, matrix_t imags){

    int m = N;
    int alignment = 32;
    int outer_array_size = (int)log2(N);

    for (int i =0; i < outer_array_size;i++){

        int inner_array_size = (int) pow(2,i);
        int aux_size = inner_array_size * 2;

        ALLOC_ALIGNED(reals[i], alignment, sizeof(data_t) * inner_array_size);
        ALLOC_ALIGNED(imags[i], alignment, sizeof(data_t) * inner_array_size);

        for (int j =0;j < inner_array_size;j++){
            reals[i][j] = cos(2.0*M_PI*(data_t)j/(data_t)aux_size);
            imags[i][j] = sin(2.0*M_PI*(data_t)j/(data_t)aux_size);
        }
    }
}

void init_look_up_inverse_d(int N, matrix_t reals, matrix_t imags){
    
    int m = N;
    int alignment = 32;
    int outer_array_size = (int)log2(N);

    for (int i =0; i < outer_array_size;i++){

        int inner_array_size = (int) pow(2,i);
        int aux_size = inner_array_size * 2;

        ALLOC_ALIGNED(reals[i], alignment, sizeof(data_t) * inner_array_size);
        ALLOC_ALIGNED(imags[i], alignment, sizeof(data_t) * inner_array_size);

        for (int j =0;j < inner_array_size;j++){
            reals[i][j] = cos(2.0*M_PI*(data_t)j/(data_t)aux_size);
            imags[i][j] = sin(-2.0*M_PI*(data_t)j/(data_t)aux_size);
        }
    }
    
    
}


void recursive_fft_d(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n){


    if (n == 4){

        data_t o0_real = in_reals[0];
        data_t o0_imag = in_imags[0];
        data_t o1_real = in_reals[stride];
        data_t o1_imag = in_imags[stride];
        data_t o2_real = in_reals[2*stride];
        data_t o2_imag = in_imags[2*stride];
        data_t o3_real = in_reals[3*stride];
        data_t o3_imag = in_imags[3*stride];
        
        out_reals[0] = o0_real + o2_real + o1_real + o3_real;
        out_imags[0] = o0_imag + o2_imag + o1_imag + o3_imag;

        out_reals[1] = o0_real - o2_real - o1_imag + o3_imag;
        out_imags[1] = o0_imag - o2_imag + o1_real - o3_real;

        out_reals[2] = o0_real + o2_real - o1_real - o3_real;
        out_imags[2] = o0_imag + o2_imag - o1_imag - o3_imag;

        out_reals[3] = o0_real - o2_real + o1_imag - o3_imag;
        out_imags[3] = o0_imag - o2_imag - o1_real + o3_real;


    }else {
        recursive_fft_d(in_reals, in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_fft_d(in_reals + stride, in_imags + stride, out_reals + n/2, out_imags + n/2, w_reals, w_imags, stride << 1, n >> 1);

        complex_t w, y_1_k, t, y_0_k;
        // twiddle factor outer array index
        int aux_num = log2(n)-1;

        for (int k =0; k < n/2; k+=4){

            w = LOAD(&w_reals[aux_num][k],&w_imags[aux_num][k]);
            y_1_k = LOAD(&out_reals[k+n/2],&out_imags[k+n/2]);
            y_0_k = LOAD(&out_reals[k],&out_imags[k]);
            t = MUL(w,y_1_k);
            
            STORE(&out_reals[k],&out_imags[k],ADD(y_0_k,t));
            STORE(&out_reals[k+n/2],&out_imags[k+n/2],SUB(y_0_k,t));

        }
    }

}

void recursive_inverse_fft_d(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n){
    
    if (n == 4){
        
        
        data_t o0_real = in_reals[0];
        data_t o0_imag = in_imags[0];
        data_t o1_real = in_reals[stride];
        data_t o1_imag = in_imags[stride];
        data_t o2_real = in_reals[2*stride];
        data_t o2_imag = in_imags[2*stride];
        data_t o3_real = in_reals[3*stride];
        data_t o3_imag = in_imags[3*stride];
        

        out_reals[0] = o0_real + o2_real + o1_real + o3_real;
        out_imags[0] = o0_imag + o2_imag + o1_imag + o3_imag;

        out_reals[1] = o0_real - o2_real + o1_imag - o3_imag;
        out_imags[1] = o0_imag - o2_imag - o1_real + o3_real;

        out_reals[2] = o0_real + o2_real - o1_real - o3_real;
        out_imags[2] = o0_imag + o2_imag - o1_imag - o3_imag;

        out_reals[3] = o0_real - o2_real - o1_imag + o3_imag;
        out_imags[3] = o0_imag - o2_imag + o1_real - o3_real;
        
        
    }else{
        recursive_inverse_fft_d(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_inverse_fft_d(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, w_reals, w_imags, stride << 1, n >> 1);

        
        complex_t w, y_1_k, t, y_0_k;
        int k;
        // twiddle factor outer array index
        int aux_num = log2(n)-1;

        for (k =0; k < n/2; k+=4){
            
            w = LOAD(&w_reals[aux_num][k],&w_imags[aux_num][k]);
            y_1_k = LOAD(&out_reals[k+n/2],&out_imags[k + n/2]);
            y_0_k = LOAD(&out_reals[k],&out_imags[k]);
            t = MUL(w,y_1_k);
            
            STORE(&out_reals[k],&out_imags[k],ADD(y_0_k,t));
            STORE(&out_reals[k+n/2],&out_imags[k+n/2],SUB(y_0_k,t));

        }

        
    }

    
}

void recursive_rfft_half_zero_d(array_t in_reals, array_t in_imags, array_t out_reals, array_t out_imags, matrix_t w_reals, matrix_t w_imags , int stride, int n){


    if (n == 4){

        data_t o0_real = in_reals[0];
        data_t o0_imag = in_imags[0];
        data_t o1_real = in_reals[stride];
        data_t o1_imag = in_imags[stride];
        
        out_reals[0] = o0_real + o1_real;
        out_imags[0] = o0_imag + o1_imag;

        out_reals[1] = o0_real - o1_imag;
        out_imags[1] = o0_imag + o1_real;

        out_reals[2] = o0_real - o1_real;
        out_imags[2] = o0_imag - o1_imag;

        out_reals[3] = o0_real + o1_imag;
        out_imags[3] = o0_imag - o1_real;


    }else {
        recursive_rfft_half_zero_d(in_reals, in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_rfft_half_zero_d(in_reals + stride, in_imags + stride, out_reals + n/2, out_imags + n/2, w_reals, w_imags, stride << 1, n >> 1);

        complex_t w, y_1_k, t, y_0_k;
        // twiddle factor outer array index
        int aux_num = log2(n)-1;
        
        
        for (int k =0; k < n/2; k+=4){

            w = LOAD(&w_reals[aux_num][k],&w_imags[aux_num][k]);
            y_1_k = LOAD(&out_reals[k+n/2],&out_imags[k+n/2]);
            y_0_k = LOAD(&out_reals[k],&out_imags[k]);
            t = MUL(w,y_1_k);
            
            STORE(&out_reals[k],&out_imags[k],ADD(y_0_k,t));
            STORE(&out_reals[k+n/2],&out_imags[k+n/2],SUB(y_0_k,t));

        }
    }

}


#undef data_t
#undef array_t
#undef matrix_t