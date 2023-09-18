//
// Created by ishan on 9/7/23.
//
#include "dichotomus.h"
#include <stdio.h>

#include <time.h>
#include <cblas.h>
#include <lapacke.h>
#include "utils.h"
#include "panopticon.h"

int main(){

    GLOBAL_TIMER(MILLISECONDS,CLOCK_MONOTONIC_RAW)

    srand(time(NULL));
    goto_set_num_threads( 16);
    openblas_set_num_threads( 16);
    omp_set_num_threads(10);
    struct matrix * A = alloc_matrix(70000,6);

    generate_simulated_data(A,-1,1,0,1);

    struct dichotomus_model * model =  dichotomus_model_create_alloc(A);

    //print_matrix(A);

    //print_vector(model->person_ability);
    //print_vector(model->item_difficulty);


    for (int i =0; i< 100; i++) {
        double l = calculate_likelihood(model);
        update_person_parameters(model);
        //standardize_vector_Zscore(model->person_ability);
        update_item_parameters(model);
        //standardize_vector_Zscore(model->item_difficulty);
        
        printf("%d th iteration: %lf \n",i,l);
    }

    //print_vector(model->person_ability);
    
    print_vector(model->item_difficulty);

    //print_vector(model->person_v_total_scores);

    //print_vector(model->item_i_total_scores);
    
    free_matrix(A);



}

// gcc test.c matrix.c -O3 -march=native -mtune=native -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib/ -lopenblas -flto -lpthread -fstrict-aliasing -funroll-loops