#include <dft.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vecmath.h>
#include <sse_math_func.h>
#include <panopticon.h> 
#include <omp.h>
#define USE_SSE2

#define W_REAL_0_8 1.0
#define W_REAL_1_8 0.70710678
#define W_REAL_2_8 0.0
#define W_REAL_3_8 -0.70710678

#define W_IMAG_0_8 0
#define W_IMAG_1_8 0.70710678
#define W_IMAG_2_8 1.0
#define W_IMAG_3_8 0.70710678

#define W_INV_REAL_0_8 1.0
#define W_INV_REAL_1_8 0.70710678
#define W_INV_REAL_2_8 0.0
#define W_INV_REAL_3_8 -0.70710678

#define W_INV_IMAG_0_8 0
#define W_INV_IMAG_1_8 -0.70710678
#define W_INV_IMAG_2_8 -1.0
#define W_INV_IMAG_3_8 -0.70710678

void init_look_up_table(int N, float ** restrict reals, float ** restrict imags){

    int m = N;
    int alignment = 32;
    int outer_array_size = (int)log2(N);

    for (int i =0; i < outer_array_size;i++){

        int inner_array_size = (int) pow(2,i);
        int aux_size = inner_array_size * 2;

        posix_memalign(&reals[i], alignment, sizeof(float) * inner_array_size);
        posix_memalign(&imags[i], alignment, sizeof(float) * inner_array_size);

        for (int j =0;j < inner_array_size;j++){
            reals[i][j] = cos(2.0*M_PI*(float)j/(float)aux_size);
            imags[i][j] = sin(2.0*M_PI*(float)j/(float)aux_size);
        }
    }
}

void init_look_up_inverse(int N, float ** restrict  reals, float ** restrict  imags){
    
    int m = N;
    int alignment = 32;
    int outer_array_size = (int)log2(N);

    for (int i =0; i < outer_array_size;i++){

        int inner_array_size = (int) pow(2,i);
        int aux_size = inner_array_size * 2;

        posix_memalign(&reals[i], alignment, sizeof(float) * inner_array_size);
        posix_memalign(&imags[i], alignment, sizeof(float) * inner_array_size);

        for (int j =0;j < inner_array_size;j++){
            reals[i][j] = cos(2.0*M_PI*(float)j/(float)aux_size);
            imags[i][j] = sin(-2.0*M_PI*(float)j/(float)aux_size);
        }
    }
    
}

void init_look_up_table_d(int N, double ** restrict reals, double ** restrict imags){

    int m = N;
    int alignment = 32;
    int outer_array_size = (int)log2(N);

    for (int i =0; i < outer_array_size;i++){

        int inner_array_size = (int) pow(2,i);
        int aux_size = inner_array_size * 2;

        posix_memalign(&reals[i], alignment, sizeof(double) * inner_array_size);
        posix_memalign(&imags[i], alignment, sizeof(double) * inner_array_size);

        for (int j =0;j < inner_array_size;j++){
            reals[i][j] = cos(2.0*M_PI*(double)j/(double)aux_size);
            imags[i][j] = sin(2.0*M_PI*(double)j/(double)aux_size);
        }
    }
}

void init_look_up_inverse_d(int N, double ** restrict  reals, double ** restrict  imags){
    
    int m = N;
    int alignment = 32;
    int outer_array_size = (int)log2(N);

    for (int i =0; i < outer_array_size;i++){

        int inner_array_size = (int) pow(2,i);
        int aux_size = inner_array_size * 2;

        posix_memalign(&reals[i], alignment, sizeof(double) * inner_array_size);
        posix_memalign(&imags[i], alignment, sizeof(double) * inner_array_size);

        for (int j =0;j < inner_array_size;j++){
            reals[i][j] = cos(2.0*M_PI*(double)j/(double)aux_size);
            imags[i][j] = sin(-2.0*M_PI*(double)j/(double)aux_size);
        }
    }
    
    
}


void recursive_fft(float * restrict in_reals,  float *  restrict in_imags,  float * restrict out_reals,  float * restrict out_imags, float ** restrict w_reals,float ** restrict w_imags , int stride, int n){


    if (n == 4){

        float o_real_0  =    in_reals[0];
        float o_imag_0  =    in_imags[0];
        float o_real_1  =    in_reals[stride];
        float o_imag_1  =    in_imags[stride];
        float o_real_2  =    in_reals[2*stride];
        float o_imag_2  =    in_imags[2*stride];
        float o_real_3  =    in_reals[3*stride];
        float o_imag_3  =    in_imags[3*stride];
        
        out_reals[0]    =    o_real_0 + o_real_2 + o_real_1 + o_real_3;
        out_imags[0]    =    o_imag_0 + o_imag_2 + o_imag_1 + o_imag_3;

        out_reals[1]    =    o_real_0 - o_real_2 - o_imag_1 + o_imag_3;
        out_imags[1]    =    o_imag_0 - o_imag_2 + o_real_1 - o_real_3;

        out_reals[2]    =    o_real_0 + o_real_2 - o_real_1 - o_real_3;
        out_imags[2]    =    o_imag_0 + o_imag_2 - o_imag_1 - o_imag_3;

        out_reals[3]    =    o_real_0 - o_real_2 + o_imag_1 - o_imag_3;
        out_imags[3]    =    o_imag_0 - o_imag_2 - o_real_1 + o_real_3;


    } else if ( n ==8){
        recursive_fft(in_reals, in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, 4);
        recursive_fft(in_reals + stride, in_imags + stride, out_reals +4, out_imags+4, w_reals, w_imags, stride << 1, 4);
        
        float o_real_0  =   out_reals[0];
        float o_imag_0  =   out_imags[0];

        float o_real_4  =   out_reals[4];
        float o_imag_4  =   out_imags[4];

        float o_real_1  =   out_reals[1];
        float o_imag_1  =   out_imags[1];

        float o_real_5  =   out_reals[5];
        float o_imag_5  =   out_imags[5];

        float o_real_2  =   out_reals[2];
        float o_imag_2  =   out_imags[2];

        float o_real_6  =   out_reals[6];
        float o_imag_6  =   out_imags[6];

        float o_real_3  =   out_reals[3];
        float o_imag_3  =   out_imags[3];

        float o_real_7  =   out_reals[7];
        float o_imag_7  =   out_imags[7];

        out_reals[0]    =   o_real_0 + o_real_4;
        out_imags[0]    =   o_imag_0 + o_imag_4;

        out_reals[4]    =   o_real_0 - o_real_4;
        out_imags[4]    =   o_imag_0 - o_imag_4;

        out_reals[1]    =   o_real_1 + W_REAL_1_8 * o_real_5 - W_IMAG_1_8 * o_imag_5;
        out_imags[1]    =   o_imag_1 + W_REAL_1_8 * o_imag_5 + W_IMAG_1_8 * o_real_5;

        out_reals[5]    =   o_real_1 - W_REAL_1_8 * o_real_5 + W_IMAG_1_8 * o_imag_5;
        out_imags[5]    =   o_imag_1 - W_REAL_1_8 * o_imag_5 - W_IMAG_1_8 * o_real_5;

        out_reals[2]    =   o_real_2 - o_imag_6;
        out_imags[2]    =   o_imag_2 + o_real_6;

        out_reals[6]    =   o_real_2 + o_imag_6;
        out_imags[6]    =   o_imag_2 - o_real_6;

        out_reals[3]    =   o_real_3 + W_REAL_3_8 * o_real_7 - W_IMAG_3_8 * o_imag_7;
        out_imags[3]    =   o_imag_3 + W_REAL_3_8 * o_imag_7 + W_IMAG_3_8 * o_real_7;

        out_reals[7]    =   o_real_3 - W_REAL_3_8 * o_real_7 + W_IMAG_3_8 * o_imag_7;
        out_imags[7]    =   o_imag_3 - W_REAL_3_8 * o_imag_7 - W_IMAG_3_8 * o_real_7; 


    }else{
        recursive_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, w_reals, w_imags, stride << 1, n >> 1);

        complex_8 w, y_1_k, t, y_0_k;
        int k;
        // twiddle factor outer array index
        int aux_num = log2(n)-1;
        

        for (k =0; k < n/2; k+=8){
            w = LOAD(&w_reals[aux_num][k],&w_imags[aux_num][k]);
           
            y_1_k = LOAD(&out_reals[k+n/2],&out_imags[k+n/2]);
            y_0_k = LOAD(&out_reals[k],&out_imags[k]);
            t = MUL(w,y_1_k);
            
            STORE(&out_reals[k],&out_imags[k],ADD(y_0_k,t));
            STORE(&out_reals[k+n/2],&out_imags[k+n/2],SUB(y_0_k,t));

        }
    }

}

void recursive_inverse_fft( float * restrict in_reals,  float *  restrict in_imags,  float * restrict out_reals,  float * restrict out_imags, float ** restrict w_reals, float ** restrict w_imags , int stride, int n){
    
    if (n == 4){
        
        float o0_real = in_reals[0];
        float o0_imag = in_imags[0];
        float o1_real = in_reals[stride];
        float o1_imag = in_imags[stride];
        float o2_real = in_reals[2*stride];
        float o2_imag = in_imags[2*stride];
        float o3_real = in_reals[3*stride];
        float o3_imag = in_imags[3*stride];
        

        out_reals[0] = o0_real + o2_real + o1_real + o3_real;
        out_imags[0] = o0_imag + o2_imag + o1_imag + o3_imag;

        out_reals[1] = o0_real - o2_real + o1_imag - o3_imag;
        out_imags[1] = o0_imag - o2_imag - o1_real + o3_real;

        out_reals[2] = o0_real + o2_real - o1_real - o3_real;
        out_imags[2] = o0_imag + o2_imag - o1_imag - o3_imag;

        out_reals[3] = o0_real - o2_real - o1_imag + o3_imag;
        out_imags[3] = o0_imag - o2_imag + o1_real - o3_real;
        
        
    }else if (n ==8){
        recursive_inverse_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, 4);
        recursive_inverse_fft(in_reals + stride,in_imags + stride, out_reals +4, out_imags+4, w_reals, w_imags, stride << 1, 4);

        float o_real_0  =   out_reals[0];
        float o_imag_0  =   out_imags[0];

        float o_real_4  =   out_reals[4];
        float o_imag_4  =   out_imags[4];

        float o_real_1  =   out_reals[1];
        float o_imag_1  =   out_imags[1];

        float o_real_5  =   out_reals[5];
        float o_imag_5  =   out_imags[5];

        float o_real_2  =   out_reals[2];
        float o_imag_2  =   out_imags[2];

        float o_real_6  =   out_reals[6];
        float o_imag_6  =   out_imags[6];

        float o_real_3  =   out_reals[3];
        float o_imag_3  =   out_imags[3];

        float o_real_7  =   out_reals[7];
        float o_imag_7  =   out_imags[7];

        out_reals[0]    =   o_real_0 + o_real_4;
        out_imags[0]    =   o_imag_0 + o_imag_4;

        out_reals[4]    =   o_real_0 - o_real_4;
        out_imags[4]    =   o_imag_0 - o_imag_4;

        out_reals[1]    =   o_real_1 + W_INV_REAL_1_8 * o_real_5 - W_INV_IMAG_1_8 * o_imag_5;
        out_imags[1]    =   o_imag_1 + W_INV_REAL_1_8 * o_imag_5 + W_INV_IMAG_1_8 * o_real_5;

        out_reals[5]    =   o_real_1 - W_INV_REAL_1_8 * o_real_5 + W_INV_IMAG_1_8 * o_imag_5;
        out_imags[5]    =   o_imag_1 - W_INV_REAL_1_8 * o_imag_5 - W_INV_IMAG_1_8 * o_real_5;

        out_reals[2]    =   o_real_2 + o_imag_6;
        out_imags[2]    =   o_imag_2 - o_real_6;

        out_reals[6]    =   o_real_2 - o_imag_6;
        out_imags[6]    =   o_imag_2 + o_real_6;

        out_reals[3]    =   o_real_3 + W_INV_REAL_3_8 * o_real_7 - W_INV_IMAG_3_8 * o_imag_7;
        out_imags[3]    =   o_imag_3 + W_INV_REAL_3_8 * o_imag_7 + W_INV_IMAG_3_8 * o_real_7;

        out_reals[7]    =   o_real_3 - W_INV_REAL_3_8 * o_real_7 + W_INV_IMAG_3_8 * o_imag_7;
        out_imags[7]    =   o_imag_3 - W_INV_REAL_3_8 * o_imag_7 - W_INV_IMAG_3_8 * o_real_7; 


    }else{
        recursive_inverse_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_inverse_fft(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, w_reals, w_imags, stride << 1, n >> 1);

        complex_8 w, y_1_k, t, y_0_k;
        int k;
        // twiddle factor outer array index
        int aux_num = log2(n)-1;

        
        for (k =0; k < n/2; k+=8){

            w = LOAD(&w_reals[aux_num][k],&w_imags[aux_num][k]);
           
            y_1_k = LOAD(&out_reals[k+n/2],&out_imags[k + n/2]);
            y_0_k = LOAD(&out_reals[k],&out_imags[k]);
            t = MUL(w,y_1_k);
            
            STORE(&out_reals[k],&out_imags[k],ADD(y_0_k,t));
            STORE(&out_reals[k+n/2],&out_imags[k+n/2],SUB(y_0_k,t));

        }
    }
}

void recursive_fft_d(double * restrict in_reals, double * restrict  in_imags, double * restrict out_reals, double * restrict out_imags, double ** restrict w_reals, double ** restrict w_imags , int stride, int n){


    if (n == 4){

        double o0_real = in_reals[0];
        double o0_imag = in_imags[0];
        double o1_real = in_reals[stride];
        double o1_imag = in_imags[stride];
        double o2_real = in_reals[2*stride];
        double o2_imag = in_imags[2*stride];
        double o3_real = in_reals[3*stride];
        double o3_imag = in_imags[3*stride];
        
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

        complex_4 w, y_1_k, t, y_0_k;
        // twiddle factor outer array index
        int aux_num = log2(n)-1;

        for (int k =0; k < n/2; k+=4){

            w = LOAD_4(&w_reals[aux_num][k],&w_imags[aux_num][k]);
            y_1_k = LOAD_4(&out_reals[k+n/2],&out_imags[k+n/2]);
            y_0_k = LOAD_4(&out_reals[k],&out_imags[k]);
            t = MUL_4(w,y_1_k);
            
            STORE_4(&out_reals[k],&out_imags[k],ADD_4(y_0_k,t));
            STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(y_0_k,t));

        }
    }

}

void naive_recursive_fft(complex double * restrict in, complex double * restrict out, complex double ** restrict w, int stride, int n){

    if (n == 4){
        complex double o0 = in[0];
        complex double o1 = in[stride];
        complex double o2 = in[2*stride];
        complex double o3 = in[3*stride];

        out[0] = o0 + o2 + o1 + o3;
        out[1] = o0 - o2 - I * (o1 - o3);
        out[2] = o0 + o2 - o1 - o3;
        out[3] = o0 - o2 + I * (o1 - o3);
    } else {

        naive_recursive_fft(in, out, w, stride << 1, n >> 1);
        naive_recursive_fft(in + stride, out + n/2, w, stride << 1, n >> 1);

        int aux_num = log2(n)-1;

        for (int k =0; k < n/2; k++){
            complex double w_n = w[aux_num][k];
            complex double t = w_n * out[k+n/2];
            complex double u = out[k];
            out[k] = u + t;
            out[k+n/2] = u - t;
            
        }
    }

}


void recursive_inverse_fft_d( double * restrict in_reals, double * restrict  in_imags, double * restrict out_reals, double * restrict out_imags, double ** restrict w_reals, double ** restrict w_imags ,  int stride, int n){
    
    if (n == 4){
        
        
        double o0_real = in_reals[0];
        double o0_imag = in_imags[0];
        double o1_real = in_reals[stride];
        double o1_imag = in_imags[stride];
        double o2_real = in_reals[2*stride];
        double o2_imag = in_imags[2*stride];
        double o3_real = in_reals[3*stride];
        double o3_imag = in_imags[3*stride];
        

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

        
        complex_4 w, y_1_k, t, y_0_k;
        int k;
        // twiddle factor outer array index
        int aux_num = log2(n)-1;

        for (k =0; k < n/2; k+=4){
            
            w = LOAD_4(&w_reals[aux_num][k],&w_imags[aux_num][k]);
            y_1_k = LOAD_4(&out_reals[k+n/2],&out_imags[k + n/2]);
            y_0_k = LOAD_4(&out_reals[k],&out_imags[k]);
            t = MUL_4(w,y_1_k);
            
            STORE_4(&out_reals[k],&out_imags[k],ADD_4(y_0_k,t));
            STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(y_0_k,t));

        }

        
    }

    
}




void recursive_rfft_half_zero_d(double * restrict in_reals, double * restrict  in_imags, double * restrict out_reals, double * restrict out_imags, double ** restrict w_reals, double ** restrict w_imags , int stride, int n){


    if (n == 4){

        double o0_real = in_reals[0];
        double o0_imag = in_imags[0];
        double o1_real = in_reals[stride];
        double o1_imag = in_imags[stride];
        
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

        complex_4 w, y_1_k, t, y_0_k;
        // twiddle factor outer array index
        int aux_num = log2(n)-1;
        
        
        for (int k =0; k < n/2; k+=4){

            w = LOAD_4(&w_reals[aux_num][k],&w_imags[aux_num][k]);
            y_1_k = LOAD_4(&out_reals[k+n/2],&out_imags[k+n/2]);
            y_0_k = LOAD_4(&out_reals[k],&out_imags[k]);
            t = MUL_4(w,y_1_k);
            
            STORE_4(&out_reals[k],&out_imags[k],ADD_4(y_0_k,t));
            STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(y_0_k,t));

        }
    }

}

void recursive_split_fft_d(double * restrict in_reals, double * restrict  in_imags, double * restrict out_reals, double * restrict out_imags, double ** restrict w_reals, double ** restrict w_imags , int stride, int n){
    
    if ( n == 2){
        
        out_reals[0] = in_reals[0] + in_reals[stride];
        out_imags[0] = in_imags[0] + in_imags[stride];

        out_reals[1] = in_reals[0] - in_reals[stride];
        out_imags[1] = in_imags[0] - in_imags[stride];
    } else if ( n == 4){
        
        double o0_real = in_reals[0];
        double o0_imag = in_imags[0];
        double o1_real = in_reals[stride];
        double o1_imag = in_imags[stride];
        double o2_real = in_reals[2*stride];
        double o2_imag = in_imags[2*stride];
        double o3_real = in_reals[3*stride];
        double o3_imag = in_imags[3*stride];
        
        out_reals[0] = o0_real + o2_real + o1_real + o3_real;
        out_imags[0] = o0_imag + o2_imag + o1_imag + o3_imag;

        out_reals[1] = o0_real - o2_real + o1_imag - o3_imag;
        out_imags[1] = o0_imag - o2_imag - o1_real + o3_real;

        out_reals[2] = o0_real + o2_real - o1_real - o3_real;
        out_imags[2] = o0_imag + o2_imag - o1_imag - o3_imag;

        out_reals[3] = o0_real - o2_real - o1_imag + o3_imag;
        out_imags[3] = o0_imag - o2_imag + o1_real - o3_real;

    } else if ( n < 16){
        recursive_split_fft_d(in_reals, in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_split_fft_d(in_reals + stride, in_imags + stride, out_reals + n/2, out_imags + n/2, w_reals, w_imags, stride << 2, n >> 2);
        recursive_split_fft_d(in_reals + 3*stride, in_imags + 3*stride, out_reals + 3*n/4, out_imags + 3*n/4, w_reals, w_imags, stride << 2, n >> 2);

        int aux_num = log2(n)-1;

        //printf("n/4 = %d\n",n/4);
        for (int k =0; k < n/4;k++){
            

            double Uk_real = out_reals[k];
            double Uk_imag = out_imags[k];
            // complex_4 Uk = LOAD_4(&out_reals[k],&out_imags[k]);

            double Zk_real = out_reals[k+n/2];
            double Zk_imag = out_imags[k+n/2];
            // complex_4 Zk = LOAD_4(&out_reals[k+n/2],&out_imags[k+n/2]);

            double Uk2_real = out_reals[k+n/4];
            double Uk2_imag = out_imags[k+n/4];
            // complex_4 Uk2 = LOAD_4(&out_reals[k+n/4],&out_imags[k+n/4]);

            double Zdk_real = out_reals[k+3*n/4];
            double Zdk_imag = out_imags[k+3*n/4];
            // complex_4 Zdk = LOAD_4(&out_reals[k+3*n/4],&out_imags[k+3*n/4]);

            double w1_real = w_reals[aux_num][k];
            double w1_imag = w_imags[aux_num][k];
            // complex_4 w1 = LOAD_4(&w_reals[aux_num][k],&w_imags[aux_num][k]);

            double w3_real = w_reals[aux_num][3*k];
            double w3_imag = w_imags[aux_num][3*k];
            // complex_4 w3 = LOAD_4(&w_reals[aux_num][3*k],&w_imags[aux_num][3*k]);


            double t1_w1_real = w1_real * Zk_real - w1_imag * Zk_imag;
            double t1_w3_real = w3_real * Zdk_real - w3_imag * Zdk_imag;
            double t1_w1_imag = w1_real * Zk_imag + w1_imag * Zk_real;
            double t1_w3_imag = w3_real * Zdk_imag + w3_imag * Zdk_real;
            // complex_4 t1_w1 = MUL_4(w1,Zk);
            // complex_4 t1_w3 = MUL_4(w3,Zdk);

            double t3_w1_real =  - t1_w1_imag;
            double t3_w1_imag = t1_w1_real;
            // complex_4 t3_w1 = CONJ_4(t1_w1);

            double t3_w3_real =  t1_w3_imag;
            double t3_w3_imag = - t1_w3_real;
            // complex_4 t3_w3 = CONJ_4(t1_w3);
            
            double t1_real = t1_w1_real + t1_w3_real;
            double t1_imag = t1_w1_imag + t1_w3_imag;
            // complex_4 t1 = ADD_4(t1_w1,t1_w3);

            double t3_real = t3_w1_real + t3_w3_real;
            double t3_imag = t3_w1_imag + t3_w3_imag;
            // complex_4 t3 = ADD_4(t3_w1,t3_w3);
            

            out_reals[k] = Uk_real + t1_real;
            out_imags[k] = Uk_imag + t1_imag;
            // STORE_4(&out_reals[k],&out_imags[k],ADD_4(Uk,t1));

            out_reals[k+n/2] = Uk_real - t1_real;
            out_imags[k+n/2] = Uk_imag - t1_imag;
            // STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(Uk,t1));

            out_reals[k+n/4] = Uk2_real - t3_real;
            out_imags[k+n/4] = Uk2_imag - t3_imag;
            // STORE_4(&out_reals[k+n/4],&out_imags[k+n/4],SUB_4(Uk2,t3));

            out_reals[k+3*n/4] = Uk2_real + t3_real;
            out_imags[k+3*n/4] = Uk2_imag + t3_imag;
            // STORE_4(&out_reals[k+3*n/4],&out_imags[k+3*n/4],ADD_4(Uk2,t3));
            
        }

    } else {

        recursive_split_fft_d(in_reals, in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_split_fft_d(in_reals + stride, in_imags + stride, out_reals + n/2, out_imags + n/2, w_reals, w_imags, stride << 2, n >> 2);
        recursive_split_fft_d(in_reals + 3*stride, in_imags + 3*stride, out_reals + 3*n/4, out_imags + 3*n/4, w_reals, w_imags, stride << 2, n >> 2);

        int aux_num = log2(n)-1;


        
        for (int k =0; k < n/4;k+=4){
            
            
            complex_4 Uk = LOAD_4(&out_reals[k],&out_imags[k]);


            complex_4 Zk = LOAD_4(&out_reals[k+n/2],&out_imags[k+n/2]);

            complex_4 Uk2 = LOAD_4(&out_reals[k+n/4],&out_imags[k+n/4]);

            complex_4 Zdk = LOAD_4(&out_reals[k+3*n/4],&out_imags[k+3*n/4]);

            complex_4 w1 = LOAD_4(&w_reals[aux_num][k],&w_imags[aux_num][k]);

            complex_4 w3 = LOAD_4(&w_reals[aux_num][3*k],&w_imags[aux_num][3*k]);

            complex_4 t1_w1 = MUL_4(w1,Zk);
            complex_4 t1_w3 = MUL_4(w3,Zdk);


            complex_4 t3_w1; 
            t3_w1.real = - t1_w1.imag;
            t3_w1.imag = t1_w1.real;  

            complex_4 t3_w3;
            t3_w3.real = t1_w3.imag;
            t3_w3.imag = - t1_w3.real;
            
            complex_4 t1 = ADD_4(t1_w1,t1_w3);

            complex_4 t3 = ADD_4(t3_w1,t3_w3);
            
            STORE_4(&out_reals[k],&out_imags[k],ADD_4(Uk,t1));

            STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(Uk,t1));

            STORE_4(&out_reals[k+n/4],&out_imags[k+n/4],SUB_4(Uk2,t3));

            STORE_4(&out_reals[k+3*n/4],&out_imags[k+3*n/4],ADD_4(Uk2,t3));


            
        }
    }
}