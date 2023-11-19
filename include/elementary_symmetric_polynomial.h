/*
Copyright (C) 2023 Ishan Hegde

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA */

#ifndef ELEM_SYM_POLY_ELEMENTARY_SYMMETRIC_INCLUDED
#define ELEM_SYM_POLY_ELEMENTARY_SYMMETRIC_INCLUDED

#include <stddef.h>
#include <common.h>

#define data_t mpfr_t
#define array_t data_t *
#define matrix_t data_t **

#define state_t elementary_symmetric_state_t

typedef struct state_t *state_t;


extern state_t  init_elementary_symmetric_state(int size, int precision);

extern void update_elementary_symmetric_state(state_t state, double * elements, int size);


extern double * compute_elementary_symmetric_polynomials(state_t state);

extern void free_elementary_symmetric_state(state_t state);


#undef state_t
#undef data_t
#undef array_t
#undef matrix_t

#endif // ELEM_SYM_POLY_ELEMENTARY_SYMMETRIC_INCLUDED