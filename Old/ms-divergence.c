/***************************************************************************
 *  Description:
 *      Locate major divergence in gene neighborhoods along provided
 *      inputs.
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-16  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <xtend/mem.h>
#include <biolibc/gff.h>
#include "gff-region.h"

int     compare_adjacent(int argc, char *argv[]);
void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    if ( argc < 3 )
	usage(argv);

    return compare_adjacent(argc, argv);
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
 *  2022-02-23  Jason Bacon Begin
 ***************************************************************************/

int     compare_adjacent(int argc, char *argv[])

{
    int             arg, old, new, old_count, new_count;
    bl_gff_hood_t   hoods[2];
    
    // Toggle old and new between 0 and 1, 1 and 0
    old = 0, new = 1;
    if ( (old_count = bl_gff_hood_load(&hoods[old], argv[1])) == 0 )
    {
	fprintf(stderr, "%s: Error reading %s: %s\n", argv[0],
		argv[1], strerror(errno));
	return EX_NOINPUT;
    }
    
    for (arg = 2; arg < argc; ++arg, old = (old + 1) % 2, new = (new + 1) % 2)
    {
	if ( (new_count = bl_gff_hood_load(&hoods[new], argv[arg])) == 0 )
	{
	    fprintf(stderr, "%s: Error reading %s: %s\n", argv[0],
		    argv[arg], strerror(errno));
	    return EX_NOINPUT;
	}
	if ( new_count != old_count )
	{
	    fprintf(stderr, "%s: Neighborhood sizes differ: %d %d\n",
		    argv[0], old_count, new_count);
	    return EX_DATAERR;
	}
	
	// FIXME: Use accessor macros
	printf("%-20s: %zu  %-20s: %zu  Common: %d\n",
	    hoods[old].species, hoods[old].gene_count,
	    hoods[new].species, hoods[new].gene_count,
	    bl_gff_hood_commonality(&hoods[old], &hoods[new]));
    }
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s neighborhood1.gff3 neighborhood2.gff3 [...]\n",
	    argv[0]);
    exit(EX_USAGE);
}
