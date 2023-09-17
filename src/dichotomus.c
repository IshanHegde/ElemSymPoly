#include "dichotomus.h"

#include <stdlib.h>


struct dichotomus_model * dichotomus_model_create_alloc(struct matrix * data){

    struct dichotomus_model * model = (struct dichotomus_model *) malloc(sizeof(struct dichotomus_model));

    if (model == NULL) return NULL;

    int num_cols = data->cols;
    int num_rows = data->rows;

    struct vector * item_difficulty = calloc_vector(num_cols);
    struct vector * person_ability = calloc_vector(num_rows);

    if (item_difficulty == NULL){

        if (person_ability != NULL) free_vector(person_ability);
        free(model);
        return NULL;
    }
    if (person_ability == NULL){

        if (item_difficulty != NULL) free_vector(item_difficulty);
        free(model);
        return NULL;
    } 

    double data_ij, item_j, person_i;

    for (int i =0; i< num_rows; i++){
        for (int j =0;j < num_cols; j++){

            data_ij = get_matrix_element(data, i, j);
            item_j = get_vector_element(item_difficulty, j);
            person_i = get_vector_element(person_ability, i);

            set_vector_element(item_difficulty,j, item_j - data_ij);
            set_vector_element(person_ability, i, person_i + data_ij);
        }
    }

    standardize_vector_Zscore(item_difficulty);
    standardize_vector_Zscore(person_ability);

    model->data = data;
    model->item_difficulty = item_difficulty;
    model->person_ability = person_ability;
    model->num_items = num_cols;
    model->num_persons = num_rows;

    return model;
}