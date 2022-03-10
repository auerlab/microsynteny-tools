#include <stdio.h>
#include <string.h>
#include <sysexits.h>
#include <xtend/file.h>     // "file.h" after import
#include <xtend/string.h>
#include <xtend/mem.h>
#include "alt-str.h"

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
 *  2022-03-08  Jason Bacon Begin
 ***************************************************************************/

int     xt_alt_str_read_line_malloc(alt_str_t *alt_str, FILE *stream)

{
    char    *buff;
    size_t  buff_size, len;
    int     status;
    
    //fprintf(stderr, "Reading alt_str...\n");
    buff_size = len = 0;
    status = xt_read_line_malloc(stream, &buff, &buff_size, &len);
    //fprintf(stderr, "%s\n", buff);
    if ( len > 0 )
	alt_str->count = strsplit(buff, &alt_str->strings, "|");
    //for (c = 0; c < alt_str->count; ++c)
    //    fprintf(stderr, "%zu %s\n", c, alt_str->strings[c]);
    return status;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      .B xt_alt_str_inhale_list()
 *      reads a list of strings from a file, one per line, into a pointer
 *      array.  Memory is allocated for the pointer array and for each
 *      string.
 *
 *      Memory should be freed using xt_free_strings(3) as soon as the
 *      strings are no longer needed.
 *
 *      Inhaling large amounts of data into arrays should generally be
 *      avoided in favor of more memory-efficient use-once-and-discard
 *      strategies, but may be advantageous for small lists of strings
 *      accessed repeatedly, or necessary for a few tasks such as sorting.
 *  
 *  Arguments:
 *      stream  FILE * from which strings are read, one per line
 *      list    Pointer to a char ** (poiner array), populated with strings
 *
 *  Returns:
 *      The number of strings read, XT_READ_IO_ERR on read error
 *
 *  Examples:
 *      FILE        *instream;
 *      alt_str_t   *list;
 *      ssize_t     alt_string_count;
 *
 *      alt_string_count = xt_alt_str_inhale_list(&list, instream);
 *      ...
 *      free(list);
 *
 *  See also:
 *      xt_alt_str_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-21  Jason Bacon Begin
 ***************************************************************************/

ssize_t xt_alt_str_inhale_list(alt_str_t **list, FILE *stream)

{
    size_t      list_size = 1024, c;
    alt_str_t   temp;
    
    if ( (*list = (alt_str_t *)xt_malloc(list_size, sizeof(*list))) == NULL )
    {
	fprintf(stderr, "load_strings(): Unable to allocate list.\n");
	return EX_UNAVAILABLE;
    }
    
    for (c = 0; xt_alt_str_read_line_malloc(&temp, stream) != EOF; ++c)
    {
	if ( c == list_size - 1 )
	{
	    list_size *= 2;
	    if ( (*list = (alt_str_t *)xt_realloc(*list, list_size, sizeof(*list))) == NULL )
	    {
		fprintf(stderr, "load_strings(): Unable to reallocate list.\n");
		return EX_UNAVAILABLE;
	    }
	}
	(*list)[c] = temp;
    }
    return c;
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
 *  2022-03-09  Jason Bacon Begin
 ***************************************************************************/

int     xt_alt_str_case_contains(alt_str_t *alt_str, char *str)

{
    size_t  c;
    
    for (c = 0; c < alt_str->count; ++c)
	if ( strcasecmp(alt_str->strings[c], str) == 0 )
	    return c;
    return -1;
}
