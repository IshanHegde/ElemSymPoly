//
// Created by ishan on 9/7/23.
//

#ifndef ITEM_MATRIX_H
#define ITEM_MATRIX_H


#include "block.h"




struct matrix{
    int rows;
    int cols;
    double *data;
    struct mem_block *block;
    size_t tda;
    int owner;
};



struct matrix * alloc_matrix(size_t size_1, size_t size_2);




inline double get_matrix_element(struct matrix * mat, int row, int col){
    return mat->data[mat->tda*row + col];
}


inline  void set_matrix_element(struct matrix * mat, int row, int col, double val){

    mat->data[mat->tda * row + col] = val;

}



void free_matrix(struct matrix * mat);

void print_matrix(struct matrix * mat);

void set_identity(struct matrix * mat);

void set_rand(struct matrix * mat);

void set_binary_rand(struct matrix * mat);

void matrix_multiply(const struct matrix * restrict A, const struct matrix * restrict B, const struct matrix * restrict C);

#endif //ITEM_MATRIX_H
