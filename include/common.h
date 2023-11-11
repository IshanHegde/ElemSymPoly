#ifndef ELEM_SYM_POLY_COMMON_INCLUDED
#define ELEM_SYM_POLY_COMMON_INCLUDED


#include <panopticon.h>
#include <math.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

#define PRECISION 256

#define ALLOC_ALIGNED(ptr, alignment, size) (void)posix_memalign((void **)&(ptr), alignment, (size))

#define ALLOC(size) malloc((size))

#endif //ELEM_SYM_POLY_COMMON_INCLUDED