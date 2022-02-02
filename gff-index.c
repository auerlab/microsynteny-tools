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

int     bl_gff_index_add_pos(bl_gff_index_t *gi, long file_pos,
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

int     bl_gff_index_seek_first_ge(bl_gff_index_t *gi, FILE *stream,
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


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      .B str2u64()
 *      is a super-fast hash function that converts a string of 8 or fewer
 *      characters to a 64-bit integer.  This is useful for storing lists
 *      of short strings, as it eliminates the need to use strdup(),
 *      strlcpy(), and strcmp() for processing.  Strings can be compared
 *      for equality using a straight integer comparison.  Strings of 7
 *      or fewer characters can still be accessed as a string by simply
 *      casting to char * for output, lexical comparison with strcmp(), etc.
 *      A string of 8 characters will not have a null-terminator.
 *
 *      The value returned varies depending on endianness.  Hence, hash
 *      values generated on one architecture will need to be byte swapped
 *      before comparison to values generated under a different endianness.
 *  
 *  Arguments:
 *      str     String to convert
 *
 *  Returns:
 *      uint64_t integer containing the characters in str
 *
 *  Examples:
 *      char        *s1 = "hello!", s2 = "Hello!";
 *      uint64_t    v1, v2;
 *      
 *      v1 = str2u64(s1);
 *      v2 = str2u64(s2);
 *      if ( v1 != v2 )
 *          printf("%s and %s are different.\n", (char *)&v1, (char *)&v2);
 *
 *  See also:
 *      strdup(3), strcmp(3), strlcpy(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-02  Jason Bacon Begin
 ***************************************************************************/

uint64_t    str2u64(const char *str)

{
    size_t      c;
    char        *p;
    uint64_t    v = 0;
    
    for (c = 0, p = (char *)&v; (c < sizeof(v)) && (str[c] != '\0'); ++c, ++p)
	*p = str[c];
    return v;
}
