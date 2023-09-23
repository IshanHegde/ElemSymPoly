#include <complex_block.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct complex_block * alloc_complex_block(size_t size){

    struct complex_block * block = (struct complex_block *) malloc(sizeof(struct complex_block));

    if (block == NULL){
        fprintf(stderr, "Failed to allocate memory for complex block structure.\n");
        return NULL;
    }

    //block->data = (double *) malloc(sizeof(double) * size);
    void * tmp_mem = NULL; 
    posix_memalign(&tmp_mem, sizeof(double), sizeof(double complex) * size);
    block->data = (double complex *) tmp_mem;

    if (block->data == NULL){
        fprintf(stderr, "Failed to allocate memory for data array.\n");
        free(block);
        return NULL;
    }

    return block;
}


struct complex_block * calloc_complex_block(size_t size){

    struct complex_block * block = (struct complex_block *) malloc(sizeof(struct complex_block));

    if (block == NULL){
        fprintf(stderr, "Failed to allocate memory for complex block structure.\n");
        return NULL;
    }

    void * tmp_mem = NULL; 
    posix_memalign(&tmp_mem, sizeof(double), sizeof(double complex) * size);
    memset(tmp_mem,0.0, sizeof(double complex) * size);
    block->data = (double complex *) tmp_mem;

    if (block->data == NULL){
        fprintf(stderr, "Failed to allocate memory for data array.\n");
        free(block);
        return NULL;
    }

    return block;

}

void free_complex_block(struct complex_block * block){

    free(block->data);
    free(block);
}