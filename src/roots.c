#include <roots.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

struct roots_of_unity * init_alloc_roots(int n){

    struct roots_of_unity * r_unity = (struct roots_of_unity *) malloc(sizeof(struct roots_of_unity ));

    if (r_unity == NULL){
        fprintf(stderr,"Failed to allocate memory for roots_of_unity structure.\n");
        return NULL;
    }

    struct complex_vector * vec = alloc_complex_vector(n);

    if (vec == NULL) {
        fprintf(stderr,"Failed to allocate memory for complex vector. \n");
        free(r_unity);
        return NULL;
    }

    int k;
    double complex exponent;
    
    for (k = 0; k < n;k++){
        exponent = ( 2 * M_PI * k * I) / n;
        set_complex_vector_element(vec,k,cexp(exponent));
    }

    r_unity->n = n;
    r_unity->roots = vec;

    return r_unity;
}

void free_roots_of_unity(struct roots_of_unity * roots ){
    free_complex_vector(roots->roots);
    free(roots);
}