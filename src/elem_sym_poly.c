//
// Created by ishan on 11/12/23.
//

#include <elem_sym_poly.h>
#include <elementary_symmetric_polynomial.h>
#include <stdlib.h>
double * compute_elem_sym_poly(double * elements, int N){

    elementary_symmetric_state_t  state = init_elementary_symmetric_state(N);

    update_elementary_symmetric_state(state, elements, N);

    double * poly = compute_elementary_symmetric_polynomials(state);

    free_elementary_symmetric_state(state);

    return poly;
}

void free_poly(double * poly){
    free(poly);
}