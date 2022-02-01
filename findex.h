
#ifndef _STDIO_H_
#include <stdio.h>
#endif

#define FINDEX_INIT { 0, 0, NULL }

#define FINDEX_OK               0
#define FINDEX_MALLOC_FAILED    -1

typedef struct
{
    size_t  array_size;
    size_t  position_count;
    long    *positions;     // Return type of ftell()
}   findex_t;

int findex_add_pos(findex_t *findex, FILE *stream);
