#ifndef RASCH_ROOTS_H
#define RASCH_ROOTS_H

#include <complex_vector.h>


struct roots_of_unity{

    int n;
    struct complex_vector * roots;

};

struct roots_of_unity * init_alloc_roots(int n);

void free_roots_of_unity(struct roots_of_unity * roots );

#endif // RASCH_ROOTS_H