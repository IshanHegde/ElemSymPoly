
#include "vector.h"
#include <stdio.h>
#include <stdlib.h>

extern double get_vector_element(struct vector * vec, int i);

extern void set_vector_element(struct vector * vec, int i, double value);

struct vector * alloc_vector(size_t size){

    struct vector * vec = (struct vector *) malloc(sizeof(struct vector));

    if (vec == NULL) return NULL;

    struct mem_block * block = alloc_mem_block(size);

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

    struct mem_block * block = calloc_mem_block(size);

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

    for (int i =0; i< size;i++){
        printf("%lf ", get_vector_element(vec, i));
    }
    printf("\n");
}

void free_vector(struct vector * vec){

    if (vec->owner){
        free_block(vec->block);
    }

    free(vec);
}