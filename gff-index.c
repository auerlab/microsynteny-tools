#include <stdio.h>
#include <stdlib.h>
#include <xtend/mem.h>
#include <xtend/math.h>
#include "gff-index.h"

/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *  
 *  Arguments:
 *
 *  Returns:
 *
 *  Examples:
 *
 *  Files:
 *
 *  Environment
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-01  Jason Bacon Begin
 ***************************************************************************/

int     gff_index_add_pos(gff_index_t *gi, long file_pos,
	    char *chr, uint64_t start, uint64_t end)

{
    if ( gi->count == gi->array_size )
    {
	gi->array_size += 65536;
	gi->file_pos = xt_realloc(gi->file_pos, gi->array_size, sizeof(*gi->file_pos));
	if ( gi->file_pos == NULL )
	    return BL_GFF_INDEX_MALLOC_FAILED;
	gi->start = xt_realloc(gi->start, gi->array_size, sizeof(*gi->start));
	if ( gi->start == NULL )
	    return BL_GFF_INDEX_MALLOC_FAILED;
	gi->end = xt_realloc(gi->end, gi->array_size, sizeof(*gi->end));
	if ( gi->end == NULL )
	    return BL_GFF_INDEX_MALLOC_FAILED;
	gi->chr = xt_realloc(gi->chr, gi->array_size, sizeof(*gi->chr));
	if ( gi->chr == NULL )
	    return BL_GFF_INDEX_MALLOC_FAILED;
    }
    gi->file_pos[gi->count] = file_pos;
    gi->start[gi->count] = start;
    gi->end[gi->count] = end;
    
    // Hash chr to a 64-bit integer, padding with 0s beyond the end
    gi->chr[gi->count] = str2u64(chr);
    //printf("%lx %s ", gi->chr[gi->count], (char *)&gi->chr[gi->count]);
    ++gi->count;
    return BL_GFF_INDEX_OK;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *  
 *  Arguments:
 *
 *  Returns:
 *
 *  Examples:
 *
 *  Files:
 *
 *  Environment
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-01  Jason Bacon Begin
 ***************************************************************************/

int     gff_index_seek_first_ge(gff_index_t *gi, FILE *stream,
	    char *chr, uint64_t end)

{
    size_t      c;
    uint64_t    v = str2u64(chr);

    // First find the segment containing our chromosome
    for (c = gi->count - 1; ((int64_t)c >= 0) && (gi->chr[c] != v); --c)
	;
    
    // Now find the first feature overlapping or after our position
    while ( ((int64_t)c >= 0) &&(gi->end[c] > end) )
	--c;
    
    if ( gi->end[c] < end )
	++c;
    return fseek(stream, gi->file_pos[c], SEEK_SET);
}


uint64_t    str2u64(const char *str)

{
    size_t      c;
    char        *p;
    uint64_t    v = 0;
    
    for (c = 0, p = (char *)&v; (c < sizeof(v)) && (str[c] != '\0'); ++c, ++p)
	*p = str[c];
    return v;
}
