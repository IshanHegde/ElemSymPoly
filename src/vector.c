
#include <vector.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern double get_vector_element(struct vector * vec, int i);

extern void set_vector_element(struct vector * vec, int i, double value);

struct vector * alloc_vector(size_t size){

    struct vector * vec = (struct vector *) malloc(sizeof(struct vector));

    if (vec == NULL) return NULL;

    struct block * block = alloc_block(size);

    if (block == NULL){
        free(vec);
        fprintf(stderr, "Failed to allocate memory for vector. \n");
        return NULL;
    }

    vec->block = block;
    vec->size = size;
    vec->data = block->data;
    vec->owner = 1;

    return vec;
}

struct vector * calloc_vector(size_t size){

    struct vector * vec = (struct vector *) malloc(sizeof(struct vector));

    if (vec == NULL) return NULL;

    struct block * block = calloc_block(size);

    if (block == NULL){
        free(vec);
        fprintf(stderr, "Failed to allocate memory for vector. \n");
        return NULL;
    }

    vec->block = block;
    vec->size = size;
    vec->data = block->data;
    vec->owner = 1;

    return vec;
}

void print_vector(struct vector * vec){

    if (vec == NULL) {
        fprintf(stderr, "Invalid vector.\n");
        return;
    }

    int size = vec->size;
    int i;

    for (i =0; i< size;i++){
        printf("%lf ", get_vector_element(vec, i));
    }
    printf("\n");
}
void standardize_vector_Zscore(struct vector * vec){

    double sum = 0.0;
    double diff, mean, std, vec_i;
    int i;

    int vec_size = vec->size;

    for (i =0; i<vec_size; i++){
        sum += get_vector_element(vec,i);
    }

    mean = sum/vec_size;

    sum = 0.0;

    for (i = 0; i<vec_size;i++){
        diff = get_vector_element(vec,i) - mean;
        sum += diff * diff;
    }

    std = sqrt(sum / vec_size);

    for (i = 0; i<vec_size; i++){
        vec_i = get_vector_element(vec,i);
        set_vector_element(vec, i, (vec_i - mean) / (std));
    }

}

void free_vector(struct vector * vec){

    if (vec->owner){
        free_block(vec->block);
    }

    free(vec);
}