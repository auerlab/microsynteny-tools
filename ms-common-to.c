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
#include "gff-hood.h"

int     common_to(int argc, char *argv[]);
void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    if ( argc < 3 )
	usage(argv);

    return common_to(argc, argv);
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

int     common_to(int argc, char *argv[])

{
    int             arg, count, common_count;
    bl_gff_region_t h1, h2, hn, *intersect, *new_intersect;
    
    bl_gff_region_init(&h1);
    bl_gff_region_init(&h2);
    bl_gff_region_init(&hn);

    if ( (common_count = bl_gff_region_load(&h1, argv[1])) == 0 )
    {
	fprintf(stderr, "%s: Error reading %s: %s\n", argv[0],
		argv[1], strerror(errno));
	return EX_NOINPUT;
    }
    if ( (common_count = bl_gff_region_load(&h2, argv[2])) == 0 )
    {
	fprintf(stderr, "%s: Error reading %s: %s\n", argv[0],
		argv[2], strerror(errno));
	return EX_NOINPUT;
    }
    intersect = bl_gff_region_intersect(&h1, &h2);
    
    for (arg = 3; arg < argc; ++arg)
    {
	if ( (count = bl_gff_region_load(&hn, argv[arg])) == 0 )
	{
	    fprintf(stderr, "%s: Error reading %s: %s\n", argv[0],
		    argv[arg], strerror(errno));
	    return EX_NOINPUT;
	}
	
	new_intersect = bl_gff_region_intersect(intersect, &hn);
	bl_gff_region_destroy(intersect);
	free(intersect);
	intersect = new_intersect;
	
	// FIXME: Use accessor macros
	printf("%-20s: %zu  Common: %zu\n",
		hn.species, hn.count, intersect->count);
	    
    }
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s neighborhood1.gff3 neighborhood2.gff3 [...]\n",
	    argv[0]);
    exit(EX_USAGE);
}
