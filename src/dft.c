#include <dft.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vecmath.h>


void init_look_up_table(int N, float * restrict reals, float * restrict imags){


    for (int i =0;i < N/2;i+=8){
        

        reals[i] = cos(2.0f * M_PI * (float)(i)/(float)(N));
        imags[i] = sin(2.0f * M_PI * (float)(i)/(float)(N));
    }

}

void init_look_up_inverse(int N, float * restrict reals, float * restrict imags){


    for (int i =0;i < N/2;i+=8){
        

        reals[i] = cos(2.0f * M_PI * (float)(i)/(float)(N));
        imags[i] = sin(-2.0f * M_PI * (float)(i)/(float)(N));
    }

}


void recursive_fft( float * restrict in_reals,  float *  restrict in_imags,  float * restrict out_reals,  float * restrict out_imags, float * restrict w_reals,float * restrict w_imags , int stride, int n){


    if (n ==2){


        float o1_real = in_reals[0];
        float o1_imag = in_imags[0];
        float o2_real = in_reals[stride];
        float o2_imag = in_imags[stride];

        out_reals[0] = o1_real + o2_real;
        out_imags[0] = o1_imag + o2_imag;
        out_reals[1] = o1_real - o2_real;
        out_imags[1] = o1_imag - o2_imag;
        
    } else if ( n == 4){

        recursive_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals + 2, out_imags + 2, w_reals, w_imags, stride << 1, n >> 1);
        

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
        
        float w_real;
        float w_imag;

        recursive_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals + 4, out_imags + 4, w_reals, w_imags, stride << 1, n >> 1);

        for (int k =0; k < 4; k++){
            w_real = w_reals[k];
            w_imag = w_imags[k];

            float t_real = w_real*out_reals[k + 4] - w_imag*out_imags[k + 4];
            float t_imag = w_real*out_imags[k + 4] + w_imag*out_reals[k + 4];

            float temp_real = out_reals[k];
            float temp_imag = out_imags[k];
            
            out_reals[k] = temp_real + t_real;
            out_imags[k] = temp_imag + t_imag;
            out_reals[k+4] = temp_real - t_real;
            out_imags[k+4] = temp_imag - t_imag;

        }
        
        
    } else if ( n == 16){
        recursive_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals + 8, out_imags + 8, w_reals, w_imags, stride << 1, n >> 1);

        complex_8 w, y_1_k, t, y_0_k;

        w = LOAD(&w_reals[0],&w_imags[0]);
        y_1_k = LOAD(&out_reals[8],&out_imags[8]);
        y_0_k = LOAD(&out_reals[0],&out_imags[0]);
        
        t = MUL(w,y_1_k);
        
        STORE(&out_reals[0],&out_imags[0],ADD(y_0_k,t));
        STORE(&out_reals[8],&out_imags[8],SUB(y_0_k,t));

    }  else{

        recursive_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, w_reals, w_imags, stride << 1, n >> 1);

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

void recursive_inverse_fft( float * restrict in_reals,  float *  restrict in_imags,  float * restrict out_reals,  float * restrict out_imags, float * restrict w_reals,float * restrict w_imags , int stride, int n){
    
    if (n ==2){


        float o1_real = in_reals[0];
        float o1_imag = in_imags[0];
        float o2_real = in_reals[stride];
        float o2_imag = in_imags[stride];

        out_reals[0] = o1_real + o2_real;
        out_imags[0] = o1_imag + o2_imag;
        out_reals[1] = o1_real - o2_real;
        out_imags[1] = o1_imag - o2_imag;
        
    } else if ( n == 4){

        recursive_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals + 2, out_imags + 2, w_reals, w_imags, stride << 1, n >> 1);
        

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

        out_reals[1] = o1_real + o3_imag;
        out_imags[1] = o1_imag - o3_real;

        out_reals[2] = o0_real - o2_real;
        out_imags[2] = o0_imag - o2_imag;

        out_reals[3] = o1_real - o3_imag;
        out_imags[3] = o1_imag + o3_real;



    } else if ( n == 8){
        
        float w_real;
        float w_imag;

        recursive_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals + 4, out_imags + 4, w_reals, w_imags, stride << 1, n >> 1);

        for (int k =0; k < 4; k++){
            w_real = w_reals[k];
            w_imag = w_imags[k];

            float t_real = w_real*out_reals[k + 4] - w_imag*out_imags[k + 4];
            float t_imag = w_real*out_imags[k + 4] + w_imag*out_reals[k + 4];

            float temp_real = out_reals[k];
            float temp_imag = out_imags[k];
            
            out_reals[k] = temp_real + t_real;
            out_imags[k] = temp_imag + t_imag;
            out_reals[k+4] = temp_real - t_real;
            out_imags[k+4] = temp_imag - t_imag;

        }
        
        
    } else if ( n == 16){
        recursive_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals + 8, out_imags + 8, w_reals, w_imags, stride << 1, n >> 1);

        complex_8 w, y_1_k, t, y_0_k;

        w = LOAD(&w_reals[0],&w_imags[0]);
        y_1_k = LOAD(&out_reals[8],&out_imags[8]);
        y_0_k = LOAD(&out_reals[0],&out_imags[0]);
        
        t = MUL(w,y_1_k);
        
        STORE(&out_reals[0],&out_imags[0],ADD(y_0_k,t));
        STORE(&out_reals[8],&out_imags[8],SUB(y_0_k,t));

    }  else{

        recursive_fft(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags, stride << 1, n >> 1);
        recursive_fft(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, w_reals, w_imags, stride << 1, n >> 1);

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
