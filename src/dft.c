#include <dft.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vecmath.h>

#include <panopticon.h> 




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
