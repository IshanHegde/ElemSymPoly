#ifndef PYRASCH_ELEMENTARY_SYMMETRIC_INCLUDED
#define PYRASCH_ELEMENTARY_SYMMETRIC_INCLUDED

#include <stddef.h>
#include <common.h>

#define data_t mpfr_t
#define array_t data_t *
#define matrix_t data_t **

#define state_t elementary_symmetric_state_t

typedef struct state_t *state_t;


extern state_t  init_elementary_symmetric_state(int size);

extern void update_elementary_symmetric_state(state_t state, double * elements, int size);


extern double * compute_elementary_symmetric_polynomials(state_t state);

extern void free_elementary_symmetric_state(state_t state);


#undef state_t
#undef data_t
#undef array_t
#undef matrix_t

#endif // PYRASCH_ELEMENTARY_SYMMETRIC_INCLUDED