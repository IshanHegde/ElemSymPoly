//
// Created by ishan on 9/7/23.
//
#include "matrix.h"
#include <stdio.h>

#include <time.h>
#include <cblas.h>
int main(){

    goto_set_num_threads( 16);
    openblas_set_num_threads( 16);
    struct mat * A = alloc_mat(1000,1000);

    set_rand(A);

    struct mat * B = alloc_mat(1000,1000);

    set_rand(B);

    struct mat * C = alloc_mat(1000,1000);
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int i =0;i<100;i++)
    {
        matrix_multiply(A,B,C);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("time %lf \n",elapsed_time);


    //set_element(matrix,3,3,288);
    //print_mat(A);
    //printf("\n");
    //print_mat(B);
    //printf("\n");
    //print_mat(C);
    printf("This is the matrix's first value: %lf", get_element(C,0,0));
    free_mat(A);
    free_mat(B);
    free_mat(C);

    return 0;
}

// gcc test.c matrix.c -O3 -march=native -mtune=native -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib/ -lopenblas -flto -lpthread -fstrict-aliasing -funroll-loops