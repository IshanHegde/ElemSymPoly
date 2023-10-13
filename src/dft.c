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


void init_look_up_table(int N, float * restrict reals, float * restrict imags){

    int m = N/2;

    for (int i =0;i < m;i++){
            reals[i] = cos(2.0*M_PI*(float)i/(float)m);
            imags[i] = sin(2.0*M_PI*(float)i/(float)m);
        }
}

void init_look_up_inverse(int N, float * restrict  reals, float * restrict  imags){
    
    int m = N/2;

    for (int i =0;i < m;i++){
            reals[i] = cos(2.0*M_PI*(float)i/(float)m);
            imags[i] = sin(-2.0*M_PI*(float)i/(float)m);
        }
    
}

void init_look_up_table_d(int N, double * restrict reals, double * restrict imags){

    int m = N/2;

    for (int i =0;i < m;i++){
            reals[i] = cos(2.0*M_PI*(double)i/(double)m);
            imags[i] = sin(2.0*M_PI*(double)i/(double)m);
        }
}

void init_look_up_inverse_d(int N, double * restrict  reals, double * restrict  imags){
    
    int m = N/2;

    for (int i =0;i < m;i++){
            reals[i] = cos(2.0*M_PI*(double)i/(double)m);
            imags[i] = sin(-2.0*M_PI*(double)i/(double)m);
        }
    
}


void recursive_fft(float * restrict in_reals,  float *  restrict in_imags,  float * restrict out_reals,  float * restrict out_imags, float * restrict w_reals,float * restrict w_imags , int stride, int n){


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

        out_reals[1] = o0_real - o2_real - o1_imag + o3_imag;
        out_imags[1] = o0_imag - o2_imag + o1_real - o3_real;

        out_reals[2] = o0_real + o2_real - o1_real - o3_real;
        out_imags[2] = o0_imag + o2_imag - o1_imag - o3_imag;

        out_reals[3] = o0_real - o2_real + o1_imag - o3_imag;
        out_imags[3] = o0_imag - o2_imag - o1_real + o3_real;


    } else if ( n ==8){
        recursive_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, 4);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals +4, out_imags+4, w_reals, w_imags, stride << 1, 4);

        complex_8 w, y, t;

        w = LOAD(&w_reals[0],&w_imags[0]);

        y = LOAD(&out_reals[0],&out_imags[0]);

        t = MUL(w,y);

        STORE(&out_reals[0],&out_imags[0],ADD(y,t));


    }else{
        recursive_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, w_reals, w_imags, stride << 1, n >> 1);

        complex_8 w, y_1_k, t, y_0_k;
        int k;

        for (k =0; k < n/2; k+=8){
            w = LOAD(&w_reals[k],&w_imags[k]);
           
            y_1_k = LOAD(&out_reals[k+n/2],&out_imags[k+n/2]);
            y_0_k = LOAD(&out_reals[k],&out_imags[k]);
            t = MUL(w,y_1_k);
            
            STORE(&out_reals[k],&out_imags[k],ADD(y_0_k,t));
            STORE(&out_reals[k+n/2],&out_imags[k+n/2],SUB(y_0_k,t));

        }
    }

}

void recursive_inverse_fft( float * restrict in_reals,  float *  restrict in_imags,  float * restrict out_reals,  float * restrict out_imags, float * restrict w_reals,float * restrict w_imags , int stride, int n){
    
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
        recursive_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, 4);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals +4, out_imags+4, w_reals, w_imags, stride << 1, 4);

        complex_8 w, y, t;

        w = LOAD(&w_reals[0],&w_imags[0]);

        y = LOAD(&out_reals[0],&out_imags[0]);

        t = MUL(w,y);

        STORE(&out_reals[0],&out_imags[0],ADD(y,t));


    }else{
        recursive_inverse_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_inverse_fft(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, w_reals, w_imags, stride << 1, n >> 1);

        complex_8 w, y_1_k, t, y_0_k;
        int k;

        for (k =0; k < n/2; k+=8){

            w = LOAD(&w_reals[k],&w_imags[k]);
           
            y_1_k = LOAD(&out_reals[k+n/2],&out_imags[k + n/2]);
            y_0_k = LOAD(&out_reals[k],&out_imags[k]);
            t = MUL(w,y_1_k);
            
            STORE(&out_reals[k],&out_imags[k],ADD(y_0_k,t));
            STORE(&out_reals[k+n/2],&out_imags[k+n/2],SUB(y_0_k,t));

        }
    }
}

void recursive_fft_d(double * restrict in_reals, double * restrict  in_imags, double * restrict out_reals, double * restrict out_imags, double * restrict w_reals, double * restrict w_imags , int stride, int n){

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


    }else{
        recursive_fft_d(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_fft_d(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, w_reals, w_imags, stride << 1, n >> 1);

        complex_4 w, y_1_k, t, y_0_k;
        int k;

        for (k =0; k < n/2; k+=4){

            w = LOAD_4(&w_reals[k],&w_imags[k]);
        
            y_1_k = LOAD_4(&out_reals[k+n/2],&out_imags[k+n/2]);
            y_0_k = LOAD_4(&out_reals[k],&out_imags[k]);
            t = MUL_4(w,y_1_k);
            
            STORE_4(&out_reals[k],&out_imags[k],ADD_4(y_0_k,t));
            STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(y_0_k,t));

        }
        
        
    }

}

void recursive_inverse_fft_d( double * restrict in_reals, double * restrict  in_imags, double * restrict out_reals, double * restrict out_imags, double * restrict w_reals, double * restrict w_imags ,  int stride, int n){
    
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
        
        
        for (k =0; k < n/2; k+=4){

            w = LOAD_4(&w_reals[k],&w_imags[k]);
            y_1_k = LOAD_4(&out_reals[k+n/2],&out_imags[k + n/2]);
            y_0_k = LOAD_4(&out_reals[k],&out_imags[k]);
            t = MUL_4(w,y_1_k);
            
            STORE_4(&out_reals[k],&out_imags[k],ADD_4(y_0_k,t));
            STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(y_0_k,t));

        }
    }

    
}

