
#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _INTTYPES_H_
#include <inttypes.h>
#endif

#define BL_GFF_INDEX_INIT { 0, 0, NULL, NULL }

#define BL_GFF_INDEX_OK             0
#define BL_GFF_INDEX_MALLOC_FAILED  -1
#define BL_GFF_INDEX_BAD_ARG        -2

typedef struct
{
    size_t  array_size;
    size_t  count;
    long    *file_pos;  // Return type of ftell()
    uint64_t    *chr;
    uint64_t    *start;
    uint64_t    *end;
}   bl_gff_index_t;

int bl_gff_index_add_pos(bl_gff_index_t *gi, long file_pos, char *chr, uint64_t start, uint64_t end);
int bl_gff_index_seek_first_ge(bl_gff_index_t *gi, FILE *stream, char *chr, uint64_t end);
uint64_t str2u64(const char *str);
