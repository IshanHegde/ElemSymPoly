#ifndef ARENA_INCLUDED
#define ARENA_INCLUDED

#include <stddef.h>

#define T Arena_T

typedef struct T* T;

extern T arena_new(void);

extern void arena_dispose(T *ap);

extern void *arena_alloc(T arena, size_t nbytes);

extern void *arena_calloc(T arena, size_t count, size_t nbytes);
extern void arena_free(T arena);

#endif // ARENA