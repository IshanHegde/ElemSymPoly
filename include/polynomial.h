#ifndef PYRASCH_POLYNOMIAL_INCLUDED
#define PYRASCH_POLYNOMIAL_INCLUDED

#include <stddef.h>
#include <common.h>

#define data_t mpfr_t
#define array_t data_t *
#define matrix_t data_t **

#define state_t poly_mul_state_t

typedef struct state_t *state_t;


extern state_t init_polynomial_mul_state(int poly_size);

extern void update_polynomial_mul_state(state_t state, array_t A, array_t B, int poly_size);

extern array_t polynomial_multiply(state_t state);

extern void free_polynomial_mul_state(state_t state);

#undef state_t
#undef data_t
#undef array_t
#undef matrix_t

#endif // PYRASCH_POLYNOMIAL_INCLUDED

