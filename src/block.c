#include <block.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct block * alloc_block(size_t size){

    struct block * block = (struct block *) malloc(sizeof(struct block));

    if (block == NULL){
        fprintf(stderr, "Failed to allocate memory for mem_block structure.\n");
        return NULL;
    }

    //block->data = (double *) malloc(sizeof(double) * size);
    void * tmp_mem = NULL; 
    posix_memalign(&tmp_mem, sizeof(float), sizeof(double) * size);
    block->data = (double *) tmp_mem;

    if (block->data == NULL){
        fprintf(stderr, "Failed to allocate memory for data array.\n");
        free(block);
        return NULL;
    }

    return block;
}

struct block * calloc_block(size_t size){

    struct block * block = (struct block *) malloc(sizeof(struct block));

    if (block == NULL){
        fprintf(stderr, "Failed to allocate memory for mem_block structure.\n");
        return NULL;
    }

    void * tmp_mem = NULL; 
    posix_memalign(&tmp_mem, sizeof(float), sizeof(double) * size);
    memset(tmp_mem,0.0, sizeof(double) * size);
    block->data = (double *) tmp_mem;

    if (block->data == NULL){
        fprintf(stderr, "Failed to allocate memory for data array.\n");
        free(block);
        return NULL;
    }

    return block;

}

void free_block(struct block * block){

    free(block->data);
    free(block);
}

