#ifndef ITEM_MATRIX_H
#define ITEM_MATRIX_H


#include <block.h>


struct matrix{
    int rows;
    int cols;
    double *data;
    struct block *block;
    size_t tda;
    int owner;
};



struct matrix * alloc_matrix(size_t size_1, size_t size_2);

struct matrix * calloc_matrix(size_t size_1, size_t size_2);


inline __attribute__((always_inline)) double get_matrix_element(struct matrix * mat, int row, int col){
    return mat->data[mat->tda*row + col];
}


inline __attribute__((always_inline)) void set_matrix_element(struct matrix * mat, int row, int col, double val){

    mat->data[mat->tda * row + col] = val;

}



void free_matrix(struct matrix * mat);

void print_matrix(struct matrix * mat);

void set_identity(struct matrix * mat);

void set_rand(struct matrix * mat);



void matrix_multiply(const struct matrix * restrict A, const struct matrix * restrict B, const struct matrix * restrict C);

#endif //ITEM_MATRIX_H
