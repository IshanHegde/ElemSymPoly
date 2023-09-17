
#ifndef _VECTOR_H_
#define _VECTOR_H_

#include "block.h"


struct vector{

    int size;
    double *data;
    struct mem_block *block;
    int owner;
};


struct vector * alloc_vector(size_t size);

struct vector * calloc_vector(size_t size);

inline double get_vector_element(struct vector * vec, int i){
    return vec->data[i];
}


inline void set_vector_element(struct vector * vec, int i, double value){
    vec->data[i] = value;
}

void standardize_vector_Zscore( struct vector * vec);

void print_vector(struct vector * vec);

void free_vector(struct vector * vec);

#endif // _VECTOR_H_