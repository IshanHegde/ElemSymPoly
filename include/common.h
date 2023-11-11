#ifndef RASCH_COMMON_H
#define RASCH_COMMON_H


#include <panopticon.h>
#include <math.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

#define PRECISION 256

#define ALLOC_ALIGNED(ptr, alignment, size) (void)posix_memalign((void **)&(ptr), alignment, (size))

#define ALLOC(size) malloc((size))

#endif // RASCH_COMMON_H