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
void    print_region_feature_names(bl_gff_region_t *region);
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

    printf("%-20s %9s %6s   %s\n", "Species", "Neighbors", "Common", "Neighbor genes");
    for (arg = 1; (arg < argc) &&
		  (common_count = bl_gff_region_load(&r1, argv[arg])) == 0;
		  ++arg)
	printf("%-20s %9c %6c\n", BL_GFF_REGION_SPECIES(&r1), '-', '-');
    if ( arg == argc )
	return EX_OK;
    printf("%-20s %9zu %6c  ", BL_GFF_REGION_SPECIES(&r1),
	   BL_GFF_REGION_COUNT(&r1) - 1, '*');
    //print_region_feature_names(&r1);
    putchar('\n');

    for (++arg; (arg < argc) &&
		  (common_count = bl_gff_region_load(&rn, argv[arg])) == 0;
		  ++arg)
	printf("%-20s %9c %6c\n", BL_GFF_REGION_SPECIES(&rn), '-', '-');
    if ( arg == argc )
	return EX_OK;
    intersect = bl_gff_region_intersect(&r1, &rn);
    printf("%-20s %9zu %6zu  ", BL_GFF_REGION_SPECIES(&rn),
	   BL_GFF_REGION_COUNT(&rn) - 1, BL_GFF_REGION_COUNT(intersect));
    print_region_feature_names(intersect);
    putchar('\n');
    
    snprintf(intersect_file, PATH_MAX + 1, "Intersects/%s-%s-%s",
	     BL_GFF_REGION_GOI(&r1), BL_GFF_REGION_SPECIES(&r1),
	     BL_GFF_REGION_SPECIES(&rn));
    
    for (++arg; arg < argc; ++arg)
    {
	if ( (count = bl_gff_region_load(&rn, argv[arg])) == 0 )
	    printf("%-20s %9c %6c\n", BL_GFF_REGION_SPECIES(&rn), '-', '-');
	else
	{
	    new_intersect = bl_gff_region_intersect(intersect, &rn);
	    bl_gff_region_free(intersect);
	    free(intersect);
	    intersect = new_intersect;
	    
	    printf("%-20s %9zu %6zu  ", BL_GFF_REGION_SPECIES(&rn),
		    BL_GFF_REGION_COUNT(&rn) - 1, BL_GFF_REGION_COUNT(intersect));
	    print_region_feature_names(intersect);
	    putchar('\n');
	    // Update GFF filename
	    strlcat(intersect_file, "-", PATH_MAX + 1);
	    strlcat(intersect_file, BL_GFF_REGION_SPECIES(&rn), PATH_MAX + 1);
	    strlcat(intersect_file, ".gff3", PATH_MAX + 1);
	}
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
 *  2022-03-01  Jason Bacon Begin
 ***************************************************************************/

void    print_region_feature_names(bl_gff_region_t *region)

{
    int         c;
    bl_gff_t    *gff;
    
    for (c = 0; c < BL_GFF_REGION_COUNT(region); ++c)
    {
	gff = &(BL_GFF_REGION_FEATURES_AE(region, c));
	printf(" %s", BL_GFF_FEATURE_NAME(gff));
    }
}

void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s neighborhood1.gff3 neighborhood2.gff3 [...]\n",
	    argv[0]);
    exit(EX_USAGE);
}
