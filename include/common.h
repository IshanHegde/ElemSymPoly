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


#ifndef ELEM_SYM_POLY_COMMON_INCLUDED
#define ELEM_SYM_POLY_COMMON_INCLUDED


#include <gmp.h>
#include <mpfr.h>
#include "globals.h"

#define PRECISION global_precision

#define ALLOC_ALIGNED(ptr, alignment, size) (void)posix_memalign((void **)&(ptr), alignment, (size))

#define ALLOC(size) malloc((size))

#endif //ELEM_SYM_POLY_COMMON_INCLUDED