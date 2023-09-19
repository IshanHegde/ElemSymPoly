
#ifndef _BLOCK_H_
#define _BLOCK_H_

#include <stddef.h>


struct mem_block{
    size_t size;
    double * data;
};

struct mem_block * alloc_mem_block(size_t size);

struct mem_block * calloc_mem_block(size_t size);

void free_block(struct mem_block * block);

#endif // _BLOCK_H_