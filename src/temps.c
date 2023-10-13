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
#include <omp.h>
#include <string.h>

int main(){

    GLOBAL_TIMER(MICROSECONDS,CLOCK_MONOTONIC_RAW)  
    omp_set_num_threads(16);
    int N =pow(2,4);
    int alignment =32;

    double * in_reals;
    double * in_imags;
    int result = posix_memalign((void**)&in_reals, alignment, N * sizeof(double));
    int result2 = posix_memalign((void**)&in_imags, alignment, N * sizeof(double));


    double * out_reals;
    double * out_imags;

    int result3 = posix_memalign((void**)&out_reals, alignment, N * sizeof(double));
    int result4 = posix_memalign((void**)&out_imags, alignment, N * sizeof(double));

    double * in_reals2;
    double * in_imags2;

    int result5 = posix_memalign((void**)&in_reals2, alignment, N * sizeof(double));
    int result6 = posix_memalign((void**)&in_imags2, alignment, N * sizeof(double));
/*
    for (int i =0;i < N;i++){
        
        in_reals[i]= rand() % 100;
        in_imags[i]= rand() % 100;
        //printf("Val: %lf + i%lf \n",in_reals[i],in_imags[i]);

    }
*/
    memset(out_reals,0,N*sizeof(double));
    memset(out_imags,0,N*sizeof(double));
    memset(in_reals2,0,N*sizeof(double));
    memset(in_imags2,0,N*sizeof(double));
    memset(in_imags,0,N*sizeof(double));

    in_reals[0] = 1;
    in_reals[1] =2;
    in_reals[2] = 3;
    in_reals[3] = 4;
    //in_imags[3] = -923;
   in_reals[4] = 5;
    in_reals[5] =6;
    in_reals[6] = 7;
    in_reals[7] = 8;
 //in_imags[3] = -2000;
    in_reals[8] = 9;
    in_reals[9] =10;
    in_reals[10] = 11;
    in_reals[11] = 12;
    in_reals[12] = 13;
    in_reals[13] =14;
    in_reals[14] = 15;
    in_reals[15] = 16; 
     /*{9.54594, 
    -7.67462 + 11.2102 i,
     -19.799 - 0.353553 i,
      -7.17462 - 10.7102 i,
       3.18198, -7.17462 + 10.7102 i, -19.799 + 0.353553 i, -7.67462 - 11.2102 i}*/

    double ** w_reals = malloc(sizeof(double*)*(log2(N)));
    double ** w_imags = malloc(sizeof(double*)*(log2(N)));
    double ** w_reals_inverse = malloc(sizeof(double*)*(log2(N)));
    double ** w_imags_inverse = malloc(sizeof(double*)*(log2(N)));
    WATCH("copy")
    init_look_up_table_d(N,w_reals,w_imags);
    init_look_up_inverse_d(N,w_reals_inverse,w_imags_inverse);
    STOP_WATCH("copy")
    //A = reverse_copy(vec);

    
    
    WATCH("dft")
    //dft(vec,A);
    recursive_fft_d(in_reals,in_imags,out_reals,out_imags,w_reals,w_imags,0,1,N);
    recursive_inverse_fft_d(out_reals,out_imags,in_reals2,in_imags2,w_reals_inverse,w_imags_inverse,0,1,N);

    STOP_WATCH("dft")



    for (int i =0;i<N;i++){
        printf("Valw: %lf + i%lf \n",in_reals2[i]/N,in_imags2[i]/N);
    }
    
    free(in_reals);
    free(in_imags);
    free(out_reals);
    free(out_imags);
    free(w_reals);
    free(w_imags);



    return 0;

}