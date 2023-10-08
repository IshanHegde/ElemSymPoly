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

void init_look_up_table(int N, float *  reals, float *  imags){
    /*
    int i;

    __m128 two_pi_n = _mm_set1_ps( 2.0f * M_PI / (float)(N) );
    
    for ( i =0;i < N/2;i+=4){
        
        __m128  index = _mm_set_ps(i,i+1,i+2,i+3);
        __m128 x_vec = _mm_mul_ps(index,two_pi_n);
        __m128 sinx, cosx;
        sincos_ps(x_vec,&sinx,&cosx);

        _mm_store_ps(&reals[i],cosx);
        _mm_store_ps(&imags[i],sinx);

    }
    */
    int i;
    for (i =0;i < N;i++){
        printf("float val: %f\n",sin(2.0f*M_PI*(float)i/(float)N));
        printf("int val: %d, %d\n",i,N);
        reals[i] = cos(2.0*M_PI*(double)i/(double)N);
        imags[i] = sin(2.0*M_PI*(double)i/(double)N);
    }
}


void init_look_up_inverse(int N, float  ** restrict reals, float ** restrict  imags){
/*
    int i;

    __m128 two_pi_n = _mm_set1_ps( -2.0f * M_PI / (float)(N) );
    
    for ( i =0;i < N/2;i+=4){
        
        __m128  index = _mm_set_ps(i,i+1,i+2,i+3);
        __m128 x_vec = _mm_mul_ps(index,two_pi_n);
        __m128 sinx, cosx;
        sincos_ps(x_vec,&sinx,&cosx);

        _mm_store_ps(&reals[i],cosx);
        _mm_store_ps(&imags[i],sinx);

    }


*/  int n_iter = log2(N)-2;
    reals = malloc(sizeof(float*)*n_iter);
    imags = malloc(sizeof(float*)*n_iter);
    printf("n_iter: %d\n",n_iter );
    for (int i = 0; i< n_iter;i++){
        int m = N/pow(2,i);
        reals[i] = malloc(sizeof(float)*m);
        imags[i] = malloc(sizeof(float)*m);

        for (int j =0;j < m/2;j++){
            reals[i][j] = cos(-2.0*M_PI*(double)j/(double)m);
            imags[i][j] = sin(-2.0*M_PI*(double)j/(double)m);
        }
    }
    

}

void init_look_up_table_d(int N, double ** restrict reals, double ** restrict imags){

    int alignment = 32;
    int n_iter = log2(N)-1;
    //reals = malloc(sizeof(double*)*n_iter);
    //imags = malloc(sizeof(double*)*n_iter);
    
    for (int i = 0; i< n_iter;i++){
        int m = N/pow(2,i);
        
        //reals[i] = malloc(sizeof(double)*m/2);
        //imags[i] = malloc(sizeof(double)*m/2);
        posix_memalign((void**)&reals[i], alignment, m * sizeof(double));
        posix_memalign((void**)&imags[i], alignment, m * sizeof(double));
        printf("m: %d\n",m);
        for (int j =0;j < m/2;j++){
            reals[i][j] = cos(2.0*M_PI*(double)j/(double)m);
            imags[i][j] = sin(2.0*M_PI*(double)j/(double)m);
        }
    }
}

void init_look_up_inverse_d(int N, double ** restrict  reals, double ** restrict  imags){
    
    int alignment = 32;
    int n_iter = log2(N)-1;
    //reals = malloc(sizeof(double*)*n_iter);
    //imags = malloc(sizeof(double*)*n_iter);
    
    for (int i = 0; i< n_iter;i++){
        int m = N/pow(2,i);
        
        //reals[i] = malloc(sizeof(double)*m/2);
        //imags[i] = malloc(sizeof(double)*m/2);
        posix_memalign((void**)&reals[i], alignment, m * sizeof(double));
        posix_memalign((void**)&imags[i], alignment, m * sizeof(double));

        for (int j =0;j < m/2;j++){
            reals[i][j] = cos(2.0*M_PI*(double)j/(double)m);
            imags[i][j] = sin(-2.0*M_PI*(double)j/(double)m);
        }
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

void recursive_inverse_fft( float *  in_reals_inverse,  float *   in_imags_inverse,  float *  out_reals_inverse,  float *  out_imags_inverse, float *  w_reals_inverse,float *  w_imags_inverese , int stride, int n){
    
    if (n ==2){


        float o1_real = in_reals_inverse[0];
        float o1_imag = in_imags_inverse[0];
        float o2_real = in_reals_inverse[stride];
        float o2_imag = in_imags_inverse[stride];

        out_reals_inverse[0] = o1_real + o2_real;
        out_imags_inverse[0] = o1_imag + o2_imag;
        out_reals_inverse[1] = o1_real - o2_real;
        out_imags_inverse[1] = o1_imag - o2_imag;

        
        
    } else if ( n == 4){

        recursive_inverse_fft(in_reals_inverse,in_imags_inverse, out_reals_inverse, out_imags_inverse, w_reals_inverse, w_imags_inverese, stride << 1, n >> 1);
        recursive_inverse_fft(in_reals_inverse + stride,in_imags_inverse + stride, out_reals_inverse + 2, out_imags_inverse + 2, w_reals_inverse, w_imags_inverese, stride << 1, n >> 1);
        

        float o0_real = out_reals_inverse[0];
        float o0_imag = out_imags_inverse[0];
        float o1_real = out_reals_inverse[1];
        float o1_imag = out_imags_inverse[1];
        float o2_real = out_reals_inverse[2];
        float o2_imag = out_imags_inverse[2];
        float o3_real = out_reals_inverse[3];
        float o3_imag = out_imags_inverse[3];
        
        

        out_reals_inverse[0] = o0_real + o2_real;
        out_imags_inverse[0] = o0_imag + o2_imag;

        out_reals_inverse[1] = o1_real + o3_imag;
        out_imags_inverse[1] = o1_imag - o3_real;

        out_reals_inverse[2] = o0_real - o2_real;
        out_imags_inverse[2] = o0_imag - o2_imag;

        out_reals_inverse[3] = o1_real - o3_imag;
        out_imags_inverse[3] = o1_imag + o3_real;

        //out_reals[0] = o0_real + o2_real;
        //out_imags[0] = o0_imag + o2_imag;

        //out_reals[1] = o1_real - o3_imag;
        //out_imags[1] = o1_imag + o3_real;

        //out_reals[2] = o0_real - o2_real;
        //out_imags[2] = o0_imag - o2_imag;

        //out_reals[3] = o1_real + o3_imag;
        //out_imags[3] = o1_imag - o3_real;

    } else if ( n == 8){
        
        float w_real;
        float w_imag;

        recursive_inverse_fft(in_reals_inverse,in_imags_inverse, out_reals_inverse, out_imags_inverse, w_reals_inverse, w_imags_inverese, stride << 1, n >> 1);
        recursive_inverse_fft(in_reals_inverse + stride,in_imags_inverse + stride, out_reals_inverse + 4, out_imags_inverse + 4, w_reals_inverse, w_imags_inverese, stride << 1, n >> 1);

        for (int k =0; k < 4; k++){
 
            w_real = w_reals_inverse[k];
            w_imag = w_imags_inverese[k];

            float t_real = w_real*out_reals_inverse[k + 4] - w_imag*out_imags_inverse[k + 4];
            float t_imag = w_real*out_imags_inverse[k + 4] + w_imag*out_reals_inverse[k + 4];

            float temp_real = out_reals_inverse[k];
            float temp_imag = out_imags_inverse[k];
            
            out_reals_inverse[k] = temp_real + t_real;
            out_imags_inverse[k] = temp_imag + t_imag;
            out_reals_inverse[k+4] = temp_real - t_real;
            out_imags_inverse[k+4] = temp_imag - t_imag;

        }
        
        
    } else if ( n == 16){
        recursive_inverse_fft(in_reals_inverse,in_imags_inverse, out_reals_inverse, out_imags_inverse, w_reals_inverse, w_imags_inverese, stride << 1, n >> 1);
        recursive_inverse_fft(in_reals_inverse + stride,in_imags_inverse + stride, out_reals_inverse + 8, out_imags_inverse + 8, w_reals_inverse, w_imags_inverese, stride << 1, n >> 1);

        complex_8 w, y_1_k, t, y_0_k;

        w = LOAD(&w_reals_inverse[0],&w_imags_inverese[0]);
        y_1_k = LOAD(&out_reals_inverse[8],&out_imags_inverse[8]);
        y_0_k = LOAD(&out_reals_inverse[0],&out_imags_inverse[0]);
        
        t = MUL(w,y_1_k);
        
        STORE(&out_reals_inverse[0],&out_imags_inverse[0],ADD(y_0_k,t));
        STORE(&out_reals_inverse[8],&out_imags_inverse[8],SUB(y_0_k,t));

    }  else{

        recursive_inverse_fft(in_reals_inverse,in_imags_inverse, out_reals_inverse, out_imags_inverse, w_reals_inverse, w_imags_inverese, stride << 1, n >> 1);
        recursive_inverse_fft(in_reals_inverse + stride,in_imags_inverse + stride, out_reals_inverse +n/2, out_imags_inverse+n/2, w_reals_inverse, w_imags_inverese, stride << 1, n >> 1);

        complex_8 w, y_1_k, t, y_0_k;
        int k;
        
        
        for (k =0; k < n/2; k+=8){
            
            w = LOAD(&w_reals_inverse[k],&w_imags_inverese[k]);
            y_1_k = LOAD(&out_reals_inverse[k+n/2],&out_imags_inverse[k + n/2]);
            y_0_k = LOAD(&out_reals_inverse[k],&out_imags_inverse[k]);
            
            t = MUL(w,y_1_k);
            
            STORE(&out_reals_inverse[k],&out_imags_inverse[k],ADD(y_0_k,t));
            STORE(&out_reals_inverse[k+n/2],&out_imags_inverse[k+n/2],SUB(y_0_k,t));
            
        }
        
    }
}

void recursive_fft_d(double * restrict in_reals, double * restrict  in_imags, double * restrict out_reals, double * restrict out_imags, double ** restrict w_reals, double ** restrict w_imags , int log2stride, int stride, int n){

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
        recursive_fft_d(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags,log2stride+1, stride << 1, n >> 1);
        recursive_fft_d(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, w_reals, w_imags,log2stride+1, stride << 1, n >> 1);

        complex_4 w, y_1_k, t, y_0_k;
        int k;

        for (k =0; k < n/2; k+=4){

            w = LOAD_4(&w_reals[log2stride][k],&w_imags[log2stride][k]);
            y_1_k = LOAD_4(&out_reals[k+n/2],&out_imags[k + n/2]);
            y_0_k = LOAD_4(&out_reals[k],&out_imags[k]);
            t = MUL_4(w,y_1_k);
            
            STORE_4(&out_reals[k],&out_imags[k],ADD_4(y_0_k,t));
            STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(y_0_k,t));

        }
    }

}


void recursive_inverse_fft_d( double * restrict in_reals, double * restrict  in_imags, double * restrict out_reals, double * restrict out_imags, double ** restrict w_reals, double ** restrict w_imags , int log2stride, int stride, int n){
    
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
        recursive_inverse_fft_d(in_reals,in_imags, out_reals, out_imags, w_reals, w_imags,log2stride+1, stride << 1, n >> 1);
        recursive_inverse_fft_d(in_reals + stride,in_imags + stride, out_reals +n/2, out_imags+n/2, w_reals, w_imags,log2stride+1, stride << 1, n >> 1);

        complex_4 w, y_1_k, t, y_0_k;
        int k;
        
        for (k =0; k < n/2; k+=4){

            w = LOAD_4(&w_reals[log2stride][k],&w_imags[log2stride][k]);
            y_1_k = LOAD_4(&out_reals[k+n/2],&out_imags[k + n/2]);
            y_0_k = LOAD_4(&out_reals[k],&out_imags[k]);
            t = MUL_4(w,y_1_k);
            
            STORE_4(&out_reals[k],&out_imags[k],ADD_4(y_0_k,t));
            STORE_4(&out_reals[k+n/2],&out_imags[k+n/2],SUB_4(y_0_k,t));

        }
    }

    
}
