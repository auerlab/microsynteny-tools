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
#include <sys/param.h>  // PATH_MAX
#include <xtend/mem.h>
#include <xtend/file.h>
#include <biolibc/gff.h>
#include "gff-region.h"

int     common_to(int argc, char *argv[]);
void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    if ( argc < 3 )
	usage(argv);

    return common_to(argc, argv);
}


/***************************************************************************
 *  Description:
 *      
 *  History: 
 *  Date        Name        Modification
 *  2022-02-23  Jason Bacon Begin
 ***************************************************************************/

int     common_to(int argc, char *argv[])

{
    int             arg, count, common_count, c;
    bl_gff_region_t r1, rn, *intersect, *new_intersect;
    char            intersect_file[PATH_MAX + 1];
    FILE            *intersect_stream;
    
    bl_gff_region_init(&r1);
    bl_gff_region_init(&rn);

    if ( (common_count = bl_gff_region_load(&r1, argv[1])) == 0 )
    {
	fprintf(stderr, "%s: Error reading %s: %s\n", argv[0],
		argv[1], strerror(errno));
	return EX_NOINPUT;
    }
    if ( (common_count = bl_gff_region_load(&rn, argv[2])) == 0 )
    {
	fprintf(stderr, "%s: Error reading %s: %s\n", argv[0],
		argv[2], strerror(errno));
	return EX_NOINPUT;
    }
    printf("%-20s %9s %6s\n", "Species", "Neighbors", "Common");
    printf("%-20s %9zu %6zu\n", BL_GFF_REGION_SPECIES(&r1),
	   BL_GFF_REGION_COUNT(&r1) - 1, BL_GFF_REGION_COUNT(&r1) - 1);
    intersect = bl_gff_region_intersect(&r1, &rn);
    printf("%-20s %9zu %6zu\n", BL_GFF_REGION_SPECIES(&rn),
	   BL_GFF_REGION_COUNT(&rn) - 1, BL_GFF_REGION_COUNT(intersect));
    
    snprintf(intersect_file, PATH_MAX + 1, "Intersects/%s-%s-%s",
	     BL_GFF_REGION_GOI(&r1), BL_GFF_REGION_SPECIES(&r1),
	     BL_GFF_REGION_SPECIES(&rn));
    
    for (arg = 3; arg < argc; ++arg)
    {
	if ( (count = bl_gff_region_load(&rn, argv[arg])) == 0 )
	{
	    fprintf(stderr, "%s: Error reading %s: %s\n", argv[0],
		    argv[arg], strerror(errno));
	    return EX_NOINPUT;
	}
	
	new_intersect = bl_gff_region_intersect(intersect, &rn);
	bl_gff_region_free(intersect);
	free(intersect);
	intersect = new_intersect;
	
	printf("%-20s %9zu %6zu\n", BL_GFF_REGION_SPECIES(&rn),
		BL_GFF_REGION_COUNT(&rn) - 1, BL_GFF_REGION_COUNT(intersect));
	
	strlcat(intersect_file, "-", PATH_MAX + 1);
	strlcat(intersect_file, BL_GFF_REGION_SPECIES(&rn), PATH_MAX + 1);
	strlcat(intersect_file, ".gff3", PATH_MAX + 1);
    }
    
    // Write GFF for intersect
    rmkdir("Intersects", 0777);
    if ( (intersect_stream = fopen(intersect_file, "w")) == NULL )
    {
	fprintf(stderr, "ms-common-to: Could not open %s for write: %s.\n",
		intersect_file, strerror(errno));
	exit(EX_CANTCREAT);
    }
    fprintf(intersect_stream, "##gff-version 3\n");
    fprintf(intersect_stream, "## %s\n", intersect_file);
    fprintf(intersect_stream, "## This file is simply a list of genes shared by multiple species.\n");
    fprintf(intersect_stream, "## There is no chrom, start, end, etc. since they will differ across species.\n");
    for (c = 0; c < BL_GFF_REGION_COUNT(intersect); ++c)
	bl_gff_write(&BL_GFF_REGION_FEATURES_AE(intersect, c),
		     intersect_stream, BL_GFF_FIELD_ALL);
    fclose(intersect_stream);

    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s neighborhood1.gff3 neighborhood2.gff3 [...]\n",
	    argv[0]);
    exit(EX_USAGE);
}
