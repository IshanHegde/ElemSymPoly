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



int main(){

    GLOBAL_TIMER(MICROSECONDS,CLOCK_MONOTONIC_RAW)  

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


    float * in_reals2;
    float * in_imags2;

    int result9 = posix_memalign((void**)&in_reals2, alignment, N * sizeof(float));
    int result10 = posix_memalign((void**)&in_imags2, alignment, N * sizeof(float));
    //complex_array * out = malloc(sizeof( complex_array ));
    //float *  out_reals = arena_calloc(arena,N,sizeof(float));
    //float *  out_imags = arena_calloc(arena,N,sizeof(float));

    float * out_reals;
    float * out_imags;

    int result3 = posix_memalign((void**)&out_reals, alignment, N * sizeof(float));
    int result4 = posix_memalign((void**)&out_imags, alignment, N * sizeof(float));

    float * w_reals;
    float * w_imags;

    int result5 = posix_memalign((void**)&w_reals, alignment, N * sizeof(float));
    int result6 = posix_memalign((void**)&w_imags, alignment, N * sizeof(float));

    float * w_reals_inverse;
    float * w_imags_inverse;

    int result7 = posix_memalign((void**)&w_reals_inverse, alignment, N * sizeof(float));
    int result8 = posix_memalign((void**)&w_imags_inverse, alignment, N * sizeof(float));


 
    for (int i =0;i < N;i++){
        
        in_reals[i]= rand()% 100 +1;
    }
  
/*
    in_reals[0] = 1;
    in_reals[1] =-19992;
    in_reals[2] = 3;
    in_reals[3] = 827;
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
    in_reals[15] = 4;*/

 

    WATCH("copy")
    init_look_up_table(N,w_reals,w_imags);
    init_look_up_inverse(N,w_reals_inverse,w_imags_inverse);
    STOP_WATCH("copy")
    //A = reverse_copy(vec);

    
    WATCH("dft")
    //dft(vec,A);
    recursive_fft(in_reals,in_imags,out_reals,out_imags,w_reals,w_imags,1,N);
    //recursive_inverse_fft(out_reals,out_imags,in_reals2,in_imags2,w_reals_inverse,w_imags_inverse,1,N);

    STOP_WATCH("dft")



    for (int i =0;i<4;i++){
        printf("Valw: %f + i%f \n",in_reals2[i],in_imags2[i]);
    }
    //print_complex_vector(A);

    return 0;

}