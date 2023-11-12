//
// Created by ishan on 11/12/23.
//

#include <elem_sym_poly.h>
#include <elementary_symmetric_polynomial.h>

double * compute_elem_sym_poly(double * elements, int N){

    elementary_symmetric_state_t  state = init_elementary_symmetric_state(N);

    update_elementary_symmetric_state(state, elements, N);

    double * poly = compute_elementary_symmetric_polynomials(state);

    free_elementary_symmetric_state(state);

    return poly;
}