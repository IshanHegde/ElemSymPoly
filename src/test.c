//
// Created by ishan on 9/7/23.
//
#include "dichotomus.h"
#include <stdio.h>

#include <time.h>
#include <cblas.h>
#include <lapacke.h>

int main(){

    srand(time(NULL));
    goto_set_num_threads( 16);
    openblas_set_num_threads( 16);
    struct matrix * A = alloc_matrix(100,3);

    set_binary_rand(A);

    struct dichotomus_model * model =  dichotomus_model_create_alloc(A);

    print_matrix(A);

    print_vector(model->item_difficulty);

    
    free_matrix(A);



}

// gcc test.c matrix.c -O3 -march=native -mtune=native -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib/ -lopenblas -flto -lpthread -fstrict-aliasing -funroll-loops