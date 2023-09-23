#ifndef RASCH_COMPLEX_VECTOR_H
#define RASCH_COMPLEX_VECTOR_H

#include <complex_block.h>

struct complex_vector{
    int size;
    double complex *data;
    struct complex_block *block;
    int owner;
};


struct complex_vector * alloc_complex_vector(size_t size);

struct complex_vector * calloc_complex_vector(size_t size);

inline __attribute__((always_inline)) double complex get_complex_vector_element(struct complex_vector * vec, int i){
    return vec->data[i];
}


inline __attribute__((always_inline)) void set_complex_vector_element(struct complex_vector * vec, int i, double complex value){
    vec->data[i] = value;
}

void print_complex_vector(struct complex_vector * vec);

void free_complex_vector(struct complex_vector * vec);



#endif // RASCH_COMPLEX_VECTOR_H