//
// Created by ishan on 9/7/23.
//
#include "dichotomus.h"
#include <stdio.h>

#include <time.h>
#include <cblas.h>
#include <lapacke.h>

#include "utils.h"

int main(){

    srand(time(NULL));
    goto_set_num_threads( 16);
    openblas_set_num_threads( 16);
    struct matrix * A = alloc_matrix(70,10);

    generate_simulated_data(A,-2,2,0,2);

    struct dichotomus_model * model =  dichotomus_model_create_alloc(A);

    //print_matrix(A);

    //print_vector(model->person_ability);
    //print_vector(model->item_difficulty);

    for (int i =0; i< 3; i++) {
        update_person_parameters(model);
        update_item_parameters(model);
    }
    //print_vector(model->person_ability);
    print_vector(model->item_difficulty);

    //print_vector(model->person_v_total_scores);

    //print_vector(model->item_i_total_scores);
    
    free_matrix(A);



}

// gcc test.c matrix.c -O3 -march=native -mtune=native -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib/ -lopenblas -flto -lpthread -fstrict-aliasing -funroll-loops