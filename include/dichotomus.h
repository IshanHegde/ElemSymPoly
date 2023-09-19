
#ifndef _DICOTOMUS_H_

#define _DICOTOMUS_H_

#include <matrix.h>
#include <vector.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>

struct dichotomus_model{

    struct matrix * data;

    struct vector * item_difficulty;

    struct vector * person_ability;

    struct vector * item_i_total_scores;

    struct vector * person_v_total_scores;

    int num_items;
    int num_persons;

};

double calculate_likelihood(struct dichotomus_model * model);

void update_item_parameters(struct dichotomus_model * model);

void update_person_parameters(struct dichotomus_model * model);

inline double person_v_estimator(struct dichotomus_model * model, double person_v_ability_estimate){

    
    double value =0.0;
    double item_i_difficulty;
    int i;

    for ( i = 0; i< model->num_items; i++){
        item_i_difficulty = get_vector_element(model->item_difficulty, i);
        value += 1.0 / (exp(item_i_difficulty - person_v_ability_estimate ) + 1.0 ); 
    }

    return value;
}

inline double item_i_estimator(struct dichotomus_model * model, double item_i_difficulty_estimate){

    double  person_v_ability;
    double value = 0.0;
    int v;

    for ( v = 0; v < model->num_persons; v++){
        person_v_ability = get_vector_element(model->person_ability, v);
        value += 1.0 / (exp(item_i_difficulty_estimate - person_v_ability) + 1.0);
        
    }

    return value;
}

inline double person_v_estimator_differential(struct dichotomus_model * model, double person_v_ability_estimate){


    double value = 0.0;
    double item_i_difficulty, exp_value;

    //#pragma omp parallel for reduction(+:value)
    for (int i = 0; i< model->num_items; i++){
        item_i_difficulty = get_vector_element(model->item_difficulty, i);
        exp_value = exp(item_i_difficulty - person_v_ability_estimate );
        value += exp_value / (( exp_value + 1.0 ) * ( exp_value + 1.0 )); 
    }

    return value;

}

inline double item_i_estimator_differential(struct dichotomus_model * model, double item_i_difficulty_estimate){

    double value = 0.0;
    double person_v_ability, exp_value;
    //#pragma omp parallel for reduction(-:value)
    for (int v = 0; v< model->num_persons; v++){
        person_v_ability = get_vector_element(model->person_ability, v);
        exp_value = exp(item_i_difficulty_estimate - person_v_ability);
        value -= exp_value / (( exp_value + 1.0 ) * ( exp_value + 1.0 )); 
    }
    return value;
}

struct dichotomus_model * dichotomus_model_create_alloc(struct matrix * data);

double get_dichotomus_likelihood(struct dichotomus_model * model);


#endif // _DICOTOMUS_H_