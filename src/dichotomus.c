#include <dichotomus.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define data_t double
#define array_t data_t *
#define matrix_t data_t **

#define state_t dichotomus_state_t

struct state_t{

    matrix_t data;

    array_t item_difficulty;
    array_t person_ability;
    array_t item_i_total_scores;
    array_t person_v_total_scores;

    int num_items;
    int num_persons;
}

static inline data_t person_v_estimator(state_t model, data_t person_v_ability_estimate){

    data_t value =0.0;
    data_t item_i_difficulty;


    for (int i = 0; i< model->num_items; i++){
        item_i_difficulty = model->item_difficulty[i];
        value += 1.0 / (exp(item_i_difficulty - person_v_ability_estimate ) + 1.0 ); 
    }

    return value;

}


static inline data_t item_i_estimator(state_t * model, data_t item_i_difficulty_estimate){

    data_t  person_v_ability, tmp_val;
    data_t value = 0.0;

    for (int v = 0; v < model->num_persons; v++){

        person_v_ability = model->person_ability[v];
        tmp_val = exp(item_i_difficulty_estimate - person_v_ability);
        value += 1.0 / ( tmp_val+ 1.0);
    }

    return value;
}

static inline data_t person_v_estimator_differential(state_t model, data_t person_v_ability_estimate){

    data_t value = 0.0;
    data_t item_i_difficulty, exp_value;


    for (int i = 0; i< model->num_items; i++){
        item_i_difficulty = model->item_difficulty[i];
        exp_value = exp(item_i_difficulty - person_v_ability_estimate );
        value += exp_value / (( exp_value + 1.0 ) * ( exp_value + 1.0 )); 
    }

    return value;
}

static inline data_t item_i_estimator_differential(state_t model, data_t item_i_difficulty_estimate){

    data_t value = 0.0;
    data_t person_v_ability, exp_value;


    for (int v = 0; v< model->num_persons; v++){
        person_v_ability = model->person_ability[v];
        exp_value = exp(item_i_difficulty_estimate - person_v_ability);
        value -= exp_value / (( exp_value + 1.0 ) * ( exp_value + 1.0 )); 
    }

    return value;
}


void update_item_parameters(struct dichotomus_model * model){

    double item_i_estimator_derivative, item_i_difficulty_estimate, item_i_difficulty_estimate_new, item_i_estimator_value, item_i_total_score;
    int num_of_items = model->num_items;
    int i;

    for (i = 0; i < num_of_items; i++){

        item_i_difficulty_estimate = get_vector_element(model->item_difficulty, i);
        item_i_total_score = get_vector_element(model->item_i_total_scores, i);

        item_i_estimator_value = item_i_estimator(model,item_i_difficulty_estimate);
        item_i_estimator_derivative = item_i_estimator_differential(model, item_i_difficulty_estimate);

        if (item_i_estimator_derivative < 0.02) continue;

        item_i_difficulty_estimate_new = item_i_difficulty_estimate - 1.0 / item_i_estimator_derivative * (item_i_estimator_value - item_i_total_score);

        set_vector_element(model->item_difficulty, i, item_i_difficulty_estimate_new);
    }


}

void update_person_parameters(struct dichotomus_model * model){


    double person_v_estimator_derivative, person_v_ability_estimate, person_v_ability_estimate_new, person_v_estimator_value, person_v_total_score;
    int num_of_persons = model->num_persons;
    int v;
    
    for (v = 0; v < num_of_persons; v++){

        person_v_ability_estimate = get_vector_element(model->person_ability, v);
        person_v_total_score = get_vector_element(model->person_v_total_scores, v);

        person_v_estimator_value = person_v_estimator(model, person_v_ability_estimate);
        person_v_estimator_derivative = person_v_estimator_differential(model, person_v_ability_estimate);

        if (person_v_estimator_derivative < 0.02) continue;
        
        person_v_ability_estimate_new = person_v_ability_estimate - 1.0 / person_v_estimator_derivative * (person_v_estimator_value - person_v_total_score);

        set_vector_element(model->person_ability, v, person_v_ability_estimate_new);
        
    }

}

double calculate_likelihood(struct dichotomus_model * model){


    double expression_1 = 0.0;
    double expression_2 = 0.0;
    double expression_3 = 0.0;
    int i, v;

    for (v = 0; v < model->num_persons; v++){
        expression_1 += get_vector_element(model->person_v_total_scores, v) * get_vector_element(model->person_ability, v);

    }

    for (i = 0; i < model->num_items; i++){
        expression_2 += get_vector_element(model->item_i_total_scores, i) * get_vector_element(model->item_difficulty, i);
    }

    for (v = 0; v < model->num_persons; v++){
        for (i = 0; i < model->num_items; i++){
            expression_3 += log(1.0 + exp(get_vector_element(model->person_ability, v) - get_vector_element(model->item_difficulty, i) ));
        }
    }

    return expression_1 - expression_2 - expression_3;
}

struct dichotomus_model * dichotomus_model_create_alloc(struct matrix * data){

    struct dichotomus_model * model = (struct dichotomus_model *) malloc(sizeof(struct dichotomus_model));

    if (model == NULL) return NULL;

    int num_cols = data->cols;
    int num_rows = data->rows;

    struct vector * item_difficulty = calloc_vector(num_cols);
    struct vector * person_ability = calloc_vector(num_rows);
    struct vector * item_i_total_scores = calloc_vector(num_cols);
    struct vector * person_v_total_scores = calloc_vector(num_rows);


    if (item_difficulty == NULL){

        if (person_ability != NULL) free_vector(person_ability);
        if (item_i_total_scores != NULL) free_vector(item_i_total_scores);
        if (person_v_total_scores != NULL) free_vector(person_v_total_scores);
        free(model);
        return NULL;
    }
    if (person_ability == NULL){

        if (item_difficulty != NULL) free_vector(item_difficulty);
        if (item_i_total_scores != NULL) free_vector(item_i_total_scores);
        if (person_v_total_scores != NULL) free_vector(person_v_total_scores);
        free(model);
        return NULL;
    }
    if (item_i_total_scores == NULL){
        if (item_difficulty != NULL) free_vector(item_difficulty);
        if (person_ability != NULL) free_vector(person_ability);
        if (person_v_total_scores != NULL) free_vector(person_v_total_scores);
        free(model);
        return NULL;
    }
    if (person_v_total_scores == NULL){
        if (item_difficulty != NULL) free_vector(item_difficulty);
        if (person_ability != NULL) free_vector(person_ability);
        if (item_i_total_scores != NULL) free_vector(item_i_total_scores);
        free(model);
        return NULL;
    }



    double data_vi, item_i, person_v, item_i_score, person_v_score;
    int v,i;

    for (v =0; v< num_rows; v++){
        for (i =0;i < num_cols; i++){

            data_vi = get_matrix_element(data, v, i);
            item_i = get_vector_element(item_difficulty, i);
            person_v = get_vector_element(person_ability, v);
            person_v_score = get_vector_element(person_v_total_scores, v);
            item_i_score = get_vector_element(item_i_total_scores, i);

            set_vector_element(item_difficulty,i, item_i - data_vi);
            set_vector_element(person_ability, v, person_v + data_vi);
            set_vector_element(person_v_total_scores, v, person_v_score + data_vi);
            set_vector_element(item_i_total_scores, i, item_i_score + data_vi);
        }
    }

    standardize_vector_Zscore(item_difficulty);
    standardize_vector_Zscore(person_ability);

    model->data = data;
    model->item_difficulty = item_difficulty;
    model->person_ability = person_ability;
    model->item_i_total_scores = item_i_total_scores;
    model->person_v_total_scores = person_v_total_scores;
    model->num_items = num_cols;
    model->num_persons = num_rows;
    
    return model;
}