//
// Created by ishan on 9/7/23.
//

#ifndef ITEM_MATRIX_H
#define ITEM_MATRIX_H


#include <stddef.h>


struct mem_block{
    size_t size;
    double * data;
};

struct mat{
    int rows;
    int cols;
    double *data;
    struct mem_block *block;
    size_t tda;
    int owner;
};

struct mem_block * alloc_mem_block(size_t size);

struct mat * alloc_mat(size_t size_1, size_t size_2);

inline double get_element(struct mat * matrix, int row, int col){
    return matrix->data[matrix->tda*row + col];
}

inline  void set_element(struct mat * matrix, int row, int col, double val){

    matrix->data[matrix->tda * row + col] = val;

}

void free_block(struct mem_block * block);

void free_mat(struct mat * matrix);

void print_mat(struct mat * matrix);

void set_identity(struct mat * matrix);

void set_rand(struct mat * matrix);

void matrix_multiply(const struct mat* restrict A, const struct mat* restrict B, const struct mat * restrict C);

#endif //ITEM_MATRIX_H
