#ifndef RASCH_COMPLEX_BLOCK_H_
#define RASCH_COMPLEX_BLOCK_H_

#include <stddef.h>
#include <complex.h>

struct complex_block{
    size_t size;
    double complex * data;
};

struct complex_block * alloc_complex_block(size_t size);

struct complex_block * calloc_complex_block(size_t size);

void free_complex_block(struct complex_block * block);

#endif // RASCH_COMPLEX_BLOCK_H_