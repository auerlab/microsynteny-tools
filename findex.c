#include <stdio.h>
#include <xtend/mem.h>
#include "findex.h"

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

int     findex_add_pos(findex_t *findex, FILE *stream)

{
    if ( findex->position_count == findex->array_size )
    {
	findex->array_size += 1024;
	findex->positions = xt_malloc(findex->array_size, sizeof(*findex->positions));
	if ( findex->positions == NULL )
	    return FINDEX_MALLOC_FAILED;
    }
    findex->positions[++findex->position_count] = ftell(stream);
    return FINDEX_OK;
}
