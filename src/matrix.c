
#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>



void free_mat(struct mat * matrix){
    if (matrix->owner){
        free_block(matrix->block);
    }

    free(matrix);
}



struct mat * alloc_mat(size_t size_1, size_t size_2){

    struct mat * matrix = (struct mat *) malloc(sizeof(struct mat));

    if (matrix == NULL) return NULL;

    struct mem_block * block = alloc_mem_block(size_1* size_2);

    if (block == NULL){
        free(matrix);
        fprintf(stderr, "Failed to allocate memory for matrix. \n");
        return NULL;
    }

    matrix->block = block;
    matrix->rows = size_1;
    matrix->cols = size_2;
    matrix->data = block->data;
    matrix->owner = 1;
    matrix->tda = size_2;

    return matrix;
}


void print_mat(struct mat * matrix){

    if (matrix == NULL) {
        fprintf(stderr, "Invalid matrix.\n");
        return;
    }

    int rows = matrix->rows;
    int cols = matrix->cols;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%lf ", get_element(matrix, i, j));
        }
        printf("\n");
    }
}

void set_identity(struct mat * matrix){

    if (matrix == NULL) {
        fprintf(stderr, "Invalid matrix.\n");
        return;
    }

    int rows = matrix->rows;
    int cols = matrix->cols;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set_element(matrix, i, j, 1.0);
        }

    }

}

void set_rand(struct mat * matrix){
    if (matrix == NULL) {
        fprintf(stderr, "Invalid matrix.\n");
        return;
    }

    int rows = matrix->rows;
    int cols = matrix->cols;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set_element(matrix, i, j, (double)rand() / RAND_MAX);
        }

    }
}

void matrix_multiply(const struct mat* restrict A, const struct mat* restrict B, const struct mat * restrict C){


    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->rows, B->cols, A->cols, 1.0, A->data, A->tda, B->data, B->tda, 0.0, C->data, C->tda);
}