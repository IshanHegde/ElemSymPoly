#include <stdio.h>
#include <time.h>
#include <cblas.h>
#include <lapacke.h>
#include <dichotomus.h>
#include <utils.h>
#include <panopticon.h>
#include <roots.h>
#include <dft.h>
#include <arena.h>

#include <string.h>

int main(){

    GLOBAL_TIMER(MICROSECONDS,CLOCK_MONOTONIC_RAW)  

    int N =pow(2,23);
    int alignment = 32;
    //complex_array * in = malloc(sizeof( complex_array ));

    //in->reals = malloc(sizeof(float)*N);
    //in->imags = calloc(N,sizeof(float));
    //complex_array * in = arena_alloc(arena,sizeof(complex_array));

    //float * in_reals = arena_calloc(arena,N,sizeof(float));
    //float *  in_imags = arena_calloc(arena,N,sizeof(float));
    double * in_reals;
    double * in_imags;
    int result = posix_memalign((void**)&in_reals, alignment, N * sizeof(double));
    int result2 = posix_memalign((void**)&in_imags, alignment, N * sizeof(double));


    double * in_reals2;
    double * in_imags2;

    int result9 = posix_memalign((void**)&in_reals2, alignment, N * sizeof(double));
    int result10 = posix_memalign((void**)&in_imags2, alignment, N * sizeof(double));
    //complex_array * out = malloc(sizeof( complex_array ));
    //float *  out_reals = arena_calloc(arena,N,sizeof(float));
    //float *  out_imags = arena_calloc(arena,N,sizeof(float));

    double * out_reals;
    double * out_imags;

    int result3 = posix_memalign((void**)&out_reals, alignment, N * sizeof(double));
    int result4 = posix_memalign((void**)&out_imags, alignment, N * sizeof(double));

    double ** w_reals;
    double ** w_imags;

    //int result5 = posix_memalign((void**)&w_reals, alignment, N * sizeof(double));
    //int result6 = posix_memalign((void**)&w_imags, alignment, N * sizeof(double));

    w_reals =  malloc(sizeof(double *) * (int) log2(N));
    w_imags =  malloc(sizeof(double *) * (int) log2(N));

    

    double ** w_reals_inverse;
    double ** w_imags_inverse;

    w_reals_inverse =  malloc(sizeof(double *) * (int) log2(N));
    w_imags_inverse =  malloc(sizeof(double *) * (int) log2(N));

    //int result7 = posix_memalign((void**)&w_reals_inverse, alignment, N * sizeof(double));
    //int result8 = posix_memalign((void**)&w_imags_inverse, alignment, N * sizeof(double));
    memset(in_reals,0,N*sizeof(double));
    memset(in_imags,0,N*sizeof(double));
    memset(out_reals,0,N*sizeof(double));
    memset(out_imags,0,N*sizeof(double));
    //memset(w_reals,0,N*sizeof(double));
    //memset(w_imags,0,N*sizeof(double));
    //memset(w_reals_inverse,0,N*sizeof(double));
    //memset(w_imags_inverse,0,N*sizeof(double));
    memset(in_reals2,0,N*sizeof(double));
    memset(in_imags2,0,N*sizeof(double));

 

    for (int i =0;i < N;i++){
        
        in_reals[i]= i;
        
    }

/*
    in_reals[0] = -20;
    in_reals[1] =2;
    in_reals[2] = 34;
    in_reals[3] = -3;
    in_reals[4] = 1;
    in_reals[5] =2;
    in_reals[6] = 3;
    in_reals[7] = 200;
    //in_imags[3] = -2000;
    in_reals[8] = 1;
    in_reals[9] =0;
    in_reals[10] = 3;
    in_reals[11] = -21;
    in_reals[12] = 1;
    in_reals[13] =2;
    in_reals[14] = -14;
    in_reals[15] = 4;
*/ 


    WATCH("copy")
    init_look_up_table_d(N,w_reals,w_imags);
    init_look_up_inverse_d(N,w_reals_inverse,w_imags_inverse);
    STOP_WATCH("copy")
    //A = reverse_copy(vec);
    
    WATCH("dft")
    //dft(vec,A);
    recursive_fft_d(in_reals,in_imags,out_reals,out_imags,w_reals,w_imags,1,N);
    recursive_inverse_fft_d(out_reals,out_imags,in_reals2,in_imags2,w_reals_inverse,w_imags_inverse,1,N);

    STOP_WATCH("dft")



    for (int i =0;i<2;i++){
        printf("Valw: %lf + i%lf \n",in_reals2[i]/N,in_imags2[i]/N);
    }
    //print_complex_vector(A);

    return 0;

}