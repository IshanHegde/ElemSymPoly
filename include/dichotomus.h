
#ifndef PYRASCH_DICOTOMUS_INCLUDED

#define PYRASCH_DICOTOMUS_INCLUDED

#include <common.h>

#define data_t mpfr_t
#define array_t data_t *
#define matrix_t data_t **

#define state_t dichotomus_state_t

typedef struct state_t *state_t;

data_t calculate_likelihood(state_t model);

void update_item_parameters(state_t model);

void update_person_parameters(state_t model);

state_t dichotomus_model_create_alloc(matrix_t data);

double get_dichotomus_likelihood(struct dichotomus_model * model);


#endif // PYRASCH_DICOTOMUS_INCLUDED