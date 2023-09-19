#include <matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>
#include <math.h>


extern void set_matrix_element(struct matrix * mat, int row, int col, double val);

extern double get_matrix_element(struct matrix * mat, int row, int col);

void free_matrix(struct matrix * mat){
    if (mat->owner){
        free_block(mat->block);
    }

    free(mat);
}



struct matrix * alloc_matrix(size_t size_1, size_t size_2){

    struct matrix * mat = (struct matrix *) malloc(sizeof(struct matrix));

    if (mat == NULL) return NULL;

    struct mem_block * block = alloc_mem_block(size_1* size_2);

    if (block == NULL){
        free(mat);
        fprintf(stderr, "Failed to allocate memory for matrix. \n");
        return NULL;
    }

    mat->block = block;
    mat->rows = size_1;
    mat->cols = size_2;
    mat->data = block->data;
    mat->owner = 1;
    mat->tda = size_2;

    return mat;
}

struct matrix * calloc_matrix(size_t size_1, size_t size_2){

    struct matrix * mat = (struct matrix *) malloc(sizeof(struct matrix));

    if (mat == NULL) return NULL;

    struct mem_block * block = calloc_mem_block(size_1* size_2);

    if (block == NULL){
        free(mat);
        fprintf(stderr, "Failed to allocate memory for matrix. \n");
        return NULL;
    }

    mat->block = block;
    mat->rows = size_1;
    mat->cols = size_2;
    mat->data = block->data;
    mat->owner = 1;
    mat->tda = size_2;

    return mat;
}

void print_matrix(struct matrix * mat){

    if (mat == NULL) {
        fprintf(stderr, "Invalid matrix.\n");
        return;
    }

    int rows = mat->rows;
    int cols = mat->cols;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%lf ", get_matrix_element(mat, i, j));
        }
        printf("\n");
    }
}

void set_identity(struct matrix * mat){

    if (mat == NULL) {
        fprintf(stderr, "Invalid matrix.\n");
        return;
    }

    int rows = mat->rows;
    int cols = mat->cols;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set_matrix_element(mat, i, j, 1.0);
        }

    }

}

void set_rand(struct matrix * mat){
    if (mat == NULL) {
        fprintf(stderr, "Invalid matrix.\n");
        return;
    }

    int rows = mat->rows;
    int cols = mat->cols;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set_matrix_element(mat, i, j, (double)rand() / RAND_MAX);
        }

    }
}


void matrix_multiply(const struct matrix* restrict A, const struct matrix* restrict B, const struct matrix * restrict C){


    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->rows, B->cols, A->cols, 1.0, A->data, A->tda, B->data, B->tda, 0.0, C->data, C->tda);
}