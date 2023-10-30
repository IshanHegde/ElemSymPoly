#ifndef RASCH_COMMON_H
#define RASCH_COMMON_H

#include <panopticon.h>
#include <math.h>
#include <stdio.h>

#define ALLOC_ALIGNED(ptr, alignment, size) posix_memalign((void **)&(ptr), alignment, (size))

#endif // RASCH_COMMON_H