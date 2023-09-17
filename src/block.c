#include "block.h"

#include <stdio.h>
#include <stdlib.h>

struct mem_block * alloc_mem_block(size_t size){

    struct mem_block * block = (struct mem_block *) malloc(sizeof(struct mem_block));

    if (block == NULL){
        fprintf(stderr, "Failed to allocate memory for mem_block structure.\n");
        return NULL;
    }

    block->data = (double *) malloc(sizeof(double) * size);

    if (block->data == NULL){
        fprintf(stderr, "Failed to allocate memory for data array.\n");
        free(block);
        return NULL;
    }

    return block;
}

struct mem_block * calloc_mem_block(size_t size){

    struct mem_block * block = (struct mem_block *) malloc(sizeof(struct mem_block));

    if (block == NULL){
        fprintf(stderr, "Failed to allocate memory for mem_block structure.\n");
        return NULL;
    }

    block->data = (double *) calloc(size,sizeof(double));

    if (block->data == NULL){
        fprintf(stderr, "Failed to allocate memory for data array.\n");
        free(block);
        return NULL;
    }

    return block;

}

void free_block(struct mem_block * block){

    free(block->data);
    free(block);
}

