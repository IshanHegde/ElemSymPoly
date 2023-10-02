#include <arena.h>

#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <assert.h>
#include <immintrin.h>
#define T Arena_T

#define THRESHOLD 10
#define CHUNK_SIZE  100*1024

struct T {
    T prev;
    char *avail;
    char *limit;
};

union align {
    int i;
    void *p;
    float f;
    double d;
    double complex ld;
    size_t u;
};

union header{
    struct T b;
    union align a;
};

static T freechunks;
static int nfree;


T arena_new(void) {
    T arena = malloc(sizeof (*arena));
    arena->prev = NULL;
    arena->limit = arena->avail = NULL;
    return arena;

}

void * arena_alloc(T arena, size_t nbytes){

    nbytes = ((nbytes + sizeof(union align) - 1)/(sizeof(union align)))*(sizeof (union align));

    while (nbytes > arena->limit - arena->avail){
        T ptr;
        char * limit;

        if ((ptr = freechunks) != NULL) {
            freechunks = freechunks->prev;
            nfree--;
            limit = ptr->limit;
        }
        else {
            size_t m = sizeof(union header) + nbytes + CHUNK_SIZE;
            ptr = _mm_malloc(m,32);
            limit = (char *)ptr + m;
        }

        *ptr = *arena;
        arena->avail = (char *)((union header *)ptr + 1);
        arena->limit = limit;
        arena->prev = ptr;
    }

    arena->avail += nbytes;
    return arena->avail - nbytes;
}

void * arena_calloc(T arena, size_t count, size_t nbytes) {
    void * ptr;

    ptr = arena_alloc(arena, count*nbytes );
    memset(ptr, '\0', count*nbytes);
    return ptr;
}

void arena_free(T arena){

    while (arena->prev) {
        struct T tmp = *arena->prev;

        if (nfree < THRESHOLD) {
            arena->prev->prev = freechunks;
            freechunks = arena->prev;
            nfree++;
            freechunks->limit = arena->limit;
        } 
        else free(arena->prev);
        *arena = tmp;
    }
}

void arena_dispose(T *ap) {
    assert(ap && *ap);
    arena_free(*ap);
    free(*ap);
    *ap = NULL;
}

#undef CHUNK_SIZE