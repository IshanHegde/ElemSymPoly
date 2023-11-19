#ifndef ELEM_SYM_POLY_COMMON_INCLUDED
#define ELEM_SYM_POLY_COMMON_INCLUDED


#include <gmp.h>
#include <mpfr.h>
#include <globals.h>

#define PRECISION global_precision

#define ALLOC_ALIGNED(ptr, alignment, size) (void)posix_memalign((void **)&(ptr), alignment, (size))

#define ALLOC(size) malloc((size))

#endif //ELEM_SYM_POLY_COMMON_INCLUDED