#include <complex_vector.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern double complex get_complex_vector_element(struct complex_vector * vec, int i);

extern void set_complex_vector_element(struct complex_vector * vec, int i, double complex value);

struct complex_vector * alloc_complex_vector(size_t size){

    struct complex_vector * vec = (struct complex_vector *) malloc(sizeof(struct complex_vector));

    if (vec == NULL) return NULL;

    struct complex_block * block = alloc_complex_block(size);

    if (block == NULL){
        free(vec);
        fprintf(stderr, "Failed to allocate memory for complex vector. \n");
        return NULL;
    }

    vec->block = block;
    vec->size = size;
    vec->data = block->data;
    vec->owner = 1;

    return vec;
}

struct complex_vector * calloc_complex_vector(size_t size){

    struct complex_vector * vec = (struct complex_vector *) malloc(sizeof(struct complex_vector));

    if (vec == NULL) return NULL;

    struct complex_block * block = calloc_complex_block(size);

    if (block == NULL){
        free(vec);
        fprintf(stderr, "Failed to allocate memory for complex vector. \n");
        return NULL;
    }

    vec->block = block;
    vec->size = size;
    vec->data = block->data;
    vec->owner = 1;

    return vec;

}

struct complex_vector * resize_complex_vector(struct complex_vector * vec, size_t new_size){

    if (vec == NULL) return NULL;

    if (vec->size == new_size) return vec;

    struct complex_block * block = resize_complex_block(vec->block, new_size);

    if (block == NULL){
        fprintf(stderr, "Failed to allocate memory for complex vector. \n");
        return NULL;
    }

    vec->block = block;
    vec->size = new_size;
    vec->data = block->data;
    vec->owner = 1;

    return vec;
}

void print_complex_vector(struct complex_vector * vec){

    if (vec == NULL) {
        fprintf(stderr, "Invalid complex vector.\n");
        return;
    }

    int size = vec->size;
    int i;
    double real, imag;
    double complex val;

    printf("[ ");
    for (i =0; i< size-1;i++){
        val =  get_complex_vector_element(vec, i);
        real = creal(val);
        imag = cimag(val);
        printf("%lf + %lfi, \n",real,imag);
        
    }

    val =  get_complex_vector_element(vec, i);
    real = creal(val);
    imag = cimag(val);
    printf("%lf + %lfi ]\n",real,imag);
}

void free_complex_vector(struct complex_vector * vec){

    if (vec->owner){
        free_complex_block(vec->block);
    }

    free(vec);
}