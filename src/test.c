//
// Created by ishan on 9/7/23.
//
#include <stdio.h>
#include <time.h>
#include <cblas.h>
#include <lapacke.h>
#include <dichotomus.h>
#include <utils.h>
#include <panopticon.h>
#include <roots.h>
#include <dft.h>

#include <inttypes.h>
#include <fftw3.h>

int main(){

    GLOBAL_TIMER(MILLISECONDS,CLOCK_MONOTONIC_RAW)
    
    int N = 16;

    struct  vector * vec = alloc_vector(N);

    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    
    //fftw_execute(p);

    

    
    for (int i =0;i < N;i++){
        in[i] = rand()% 1000;
        set_vector_element(vec,i, rand()% 1000);
    }
    
    
    
    WATCH("fftw")

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p);
        
    STOP_WATCH("fftw")
    WATCH("dft")
    struct  complex_vector * a =  dft(vec);
    STOP_WATCH("dft")
    

    for (int i =0;i < N;i++){
        printf("A: %lf I: %lf \n",creal(out[i]),cimag(out[i]));
    }
    
    print_complex_vector(a);
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
   



    return 0;    

}
// gcc test.c matrix.c -O3 -march=native -mtune=native -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib/ -lopenblas -flto -lpthread -fstrict-aliasing -funroll-loops