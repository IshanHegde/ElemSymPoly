#include <stdio.h>
#include <time.h>
#include <cblas.h>
#include <lapacke.h>
#include <dichotomus.h>
#include <utils.h>
#include <panopticon.h>
#include <roots.h>
#include <dft.h>
#include <vecmath.h>
#include <arena.h>

#define W(N,k) (cexp(-2.0f * M_PI * I * (float)(k)/(float)(N)))


typedef struct complex_array{

    float * reals;
    float * imags;

} complex_array;

complex_array * get_even(complex_array * A,Arena_T arena, int n){

    //complex_array * ret = malloc(sizeof(complex_array));
    complex_array * ret = arena_alloc(arena,sizeof(complex_array));
    //float * reals = calloc(n/2,sizeof(float));
    //float * imags = calloc(n/2,sizeof(float));
    float * reals = arena_alloc(arena,sizeof(float)*n/2);
    float * imags = arena_alloc(arena,sizeof(float)*n/2);

    ret->reals = reals;
    ret->imags = imags;

    for (int i =0;i < n/2;i++){
        ret->reals[i] = A->reals[2*i];
        ret->imags[i] = A->imags[2*i];
    }
    

    

    return ret;
}

complex_array * get_odd(complex_array * A, Arena_T arena, int n){

    //complex_array * ret = malloc(sizeof(complex_array));
    //float * reals = calloc(n/2,sizeof(float));
    //float * imags = calloc(n/2,sizeof(float));
    complex_array * ret = arena_alloc(arena,sizeof(complex_array));
    float * reals = arena_alloc(arena,sizeof(float)*n/2);
    float * imags = arena_alloc(arena,sizeof(float)*n/2);
    ret->reals = reals;
    ret->imags = imags;

    for (int i =0;i < n/2;i++){
        ret->reals[i] = A->reals[2*i+1];
        ret->imags[i] = A->imags[2*i+1];
    }

    return ret;
}



complex_array * recursive_fft(complex_array * A, complex_array * look_up_table, Arena_T arena, int n){

    //printf("N: %d\n",n);
    if (n ==1){
        return A;
    }
    else if (n ==2){

        float o1_real = A->reals[0];
        float o1_imag = A->imags[0];
        float o2_real = A->reals[1];
        float o2_imag = A->imags[1];

        A->reals[0] = o1_real + o2_real;
        A->imags[0] = o1_imag + o2_imag;
        A->reals[1] = o1_real - o2_real;
        A->imags[1] = o1_imag - o2_imag;
        

        return A;
    } else if ( n == 4){
        
        float o0_real = A->reals[0];
        float o0_imag = A->imags[0];
        float o1_real = A->reals[1];
        float o1_imag = A->imags[1];
        float o2_real = A->reals[2];
        float o2_imag = A->imags[2];
        float o3_real = A->reals[3];
        float o3_imag = A->imags[3];


        //A[0] = o0 + o1 + o2 + o3;
        A->reals[0] = o0_real + o1_real + o2_real + o3_real;
        A->imags[0] = o0_imag + o1_imag + o2_imag + o3_imag;
        //A[1] = o0 + I*o1 - o2 - I*o3;
        A->reals[1] = o0_real - o1_imag - o2_real + o3_imag;
        A->imags[1] = o0_imag + o1_real - o2_imag - o3_real;
        //A[2] = o0 - o1 + o2  -o3;
        A->reals[2] = o0_real - o1_real + o2_real  -o3_real;
        A->imags[2] = o0_imag - o1_imag + o2_imag  -o3_imag;
        //A[3] = o0 - I*o1 - o2 + I*o3;
        A->reals[3] = o0_real + o1_imag - o2_real - o3_imag;
        A->imags[3] = o0_imag - o1_real - o2_imag + o3_real;
        
        return A;
    } else if ( n == 8){
        WATCH("8")
        float w_n_real = cos(-2.0f * M_PI * (float)(1)/(float)(n));
        float w_n_imag = sin(-2.0f * M_PI * (float)(1)/(float)(n));
        float w_real =1;
        float w_imag =0;

        complex_array * A_0 = get_even(A,arena,n);
        complex_array * A_1 = get_odd(A,arena,n);

        complex_array * y_0 = recursive_fft(A_0,look_up_table,arena, n/2);
        complex_array * y_1 = recursive_fft(A_1,look_up_table,arena, n/2);

        //complex_array * y_k = malloc(sizeof(complex_array ));
        //y_k->reals = calloc(n,sizeof(float));
        //y_k->imags = calloc(n,sizeof(float));
        complex_array * y_k = arena_alloc(arena,sizeof(complex_array));
        y_k->reals = arena_alloc(arena,sizeof(float)*n);
        y_k->imags = arena_alloc(arena,sizeof(float)*n);

        for (int k =0; k < n/2; k++){
            float t_real = w_real*y_1->reals[k] - w_imag*y_1->imags[k];
            float t_imag = w_real*y_1->imags[k] + w_imag*y_1->reals[k];

            
            y_k->reals[k] = y_0->reals[k] + t_real;
            y_k->imags[k] = y_0->imags[k] + t_imag;
            y_k->reals[k+n/2] = y_0->reals[k] - t_real;
            y_k->imags[k+n/2] = y_0->imags[k] - t_imag;


            w_real =  w_real*w_n_real - w_imag*w_n_imag;
            w_imag =  w_real*w_n_imag + w_imag*w_n_real;
        }
        STOP_WATCH("8")
        return y_k;
    } else {

        
        complex_array * A_0 = get_even(A,arena,n);
        complex_array * A_1 = get_odd(A,arena,n);
        

        complex_array * y_0 = recursive_fft(A_0,look_up_table,arena, n >> 1);
        complex_array * y_1 = recursive_fft(A_1,look_up_table,arena, n >> 1);

        
        //complex_array * y_k = malloc(sizeof(complex_array ));
        //y_k->reals = calloc(n,sizeof(float));
        //y_k->imags = calloc(n,sizeof(float));
        complex_array * y_k = arena_alloc(arena,sizeof(complex_array));
        y_k->reals = arena_alloc(arena,sizeof(float)*n);
        y_k->imags = arena_alloc(arena,sizeof(float)*n);
        
        complex_8 w, y_1_k, t, y_0_k, y_k_k, y_k_k_n2;
        int k;
        WATCH("loop")
        //#pragma omp parallel for schedule(static) num_threads(4) private(k, w, y_1_k, t, y_0_k, y_k_k, y_k_k_n2, n)
        for (k =0; k < n/2; k+=8){

             w = LOAD(&look_up_table->reals[k],&look_up_table->imags[k]);

            
             y_1_k = LOAD(&y_1->reals[k],&y_1->imags[k]);
             t = MUL(w,y_1_k);

             y_0_k = LOAD(&y_0->reals[k],&y_0->imags[k]);

             y_k_k = ADD(y_0_k,t);
            STORE(&y_k->reals[k],&y_k->imags[k],y_k_k);

             y_k_k_n2 = SUB(y_0_k,t);
            STORE(&y_k->reals[k+n/2],&y_k->imags[k+n/2],y_k_k_n2);

        }
        STOP_WATCH("loop")
        
        return y_k;
    }

}



complex_array * init_look_up_table(int N, Arena_T arena){

    //complex_array * look_up_table = malloc(sizeof(complex_array));
    //look_up_table->reals = malloc(sizeof(float)*N);
    //look_up_table->imags = malloc(sizeof(float)*N);

    complex_array * look_up_table = arena_alloc(arena,sizeof(complex_array));
    look_up_table->reals = arena_alloc(arena,sizeof(float)*N);
    look_up_table->imags = arena_alloc(arena,sizeof(float)*N);

    for (int i =0;i < N/2;i++){
        look_up_table->reals[i] = cos(-2.0f * M_PI * (float)(i)/(float)(N));
        look_up_table->imags[i] = sin(-2.0f * M_PI * (float)(i)/(float)(N));
    }

    return look_up_table;
}

int main(){

    GLOBAL_TIMER(MILLISECONDS,CLOCK_MONOTONIC_RAW)  
    Arena_T arena = arena_new();
    int N =65536;

    //complex_array * in = malloc(sizeof( complex_array ));

    //in->reals = malloc(sizeof(float)*N);
    //in->imags = calloc(N,sizeof(float));
    complex_array * in = arena_alloc(arena,sizeof(complex_array));
    in->reals = arena_alloc(arena,sizeof(float)*N);
    in->imags = arena_alloc(arena,sizeof(float)*N);

    for (int i =0;i < N;i++){
        
        in->reals[i]= rand()% 100 +1;
    }

/*
    in->reals[0] = 1;
    in->reals[1] =2;
    in->reals[2] = 3;
    in->reals[3] = -114;
*/
    

    WATCH("copy")
    complex_array * look_up_table = init_look_up_table(N,arena);
    STOP_WATCH("copy")
    //A = reverse_copy(vec);

    
    WATCH("dft")
    //dft(vec,A);
    complex_array * y = recursive_fft(in,look_up_table,arena,N);

    STOP_WATCH("dft")


    for (int i =0;i<4;i++){
        printf("Val: %f + i%f \n",y->reals[i],y->imags[i]);
    }
    //print_complex_vector(A);
    arena_dispose(&arena);
    return 0;

}