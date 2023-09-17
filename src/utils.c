#include "utils.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void generate_simulated_data(struct matrix * mat, double min_difficulty, double max_difficulty, double ability_mean, double ability_std){

    if (mat == NULL) {
        fprintf(stderr, "Invalid matrix.\n");
        return;
    }

    double threshold = 0.0;
    double mean = ability_mean; // Mean of the normal distribution
    double standard_deviation = ability_std;

    

    int rows = mat->rows;
    int cols = mat->cols;

    double difficulty_increment = (max_difficulty - min_difficulty) / (cols - 1.0);

    struct vector * difficulty = calloc_vector(cols);
    struct vector * ability = calloc_vector(rows);

    for (int i = 0; i < cols; i++){
        double diff_val =  min_difficulty + ( (double) i*difficulty_increment);
       
        set_vector_element(difficulty,i, diff_val);
    }

    for (int v = 0; v < rows; v++){
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;
        double z = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
        set_vector_element(ability,v, mean + standard_deviation * z);
    }

    //print_vector(ability);
    print_vector(difficulty);

    for (int v = 0; v < rows; v++) {
        
        for (int i = 0; i < cols; i++) {
            
            

            double score = (get_vector_element(ability,v) - get_vector_element(difficulty,i) > threshold) ? 1.0 : 0.0;

            set_matrix_element(mat,v,i,score);

        }

    }

    free_vector(ability);
    free_vector(difficulty);
}