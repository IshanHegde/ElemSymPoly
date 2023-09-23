
#ifndef _BLOCK_H_
#define _BLOCK_H_

#include <stddef.h>


struct block{
    size_t size;
    double * data;
};

struct block * alloc_block(size_t size);

struct block * calloc_block(size_t size);

void free_block(struct block * block);

#endif // _BLOCK_H_