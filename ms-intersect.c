/***************************************************************************
 *  Description:
 *      Locate major divergence in gene regions along provided
 *      inputs.  Each GFF3 file given should contain genes surrounding
 *      the same GOI for another species.
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
#include <sys/param.h>      // PATH_MAX
#include <xtend/mem.h>
#include <xtend/file.h>
#include <xtend/string.h>   // strlcat() on Linux
#include <biolibc/gff.h>
#include "gff-region.h"

int     intersect(int argc, char *argv[]);
void    print_region_feature_names(bl_gff_region_t *region);
void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    if ( argc < 3 )
	usage(argv);

    return intersect(argc, argv);
}


/***************************************************************************
 *  Description:
 *      
 *  History: 
 *  Date        Name        Modification
 *  2022-02-23  Jason Bacon Begin
 ***************************************************************************/

int     intersect(int argc, char *argv[])

{
    int             arg, count, c, old_intersect_count, new_intersect_count, max_count;
    bl_gff_region_t r1, rn, *intersect, *new_intersect, *div_intersect;
    char            intersect_file[PATH_MAX + 1];
    FILE            *intersect_stream;
    
    bl_gff_region_init(&r1);
    bl_gff_region_init(&rn);

    /*
     *  Load first two genes, skipping non-existent files
     */

    printf("%-20s %2s %2s %2s  %s\n", "Species", "Ne", "Co",
	    "Ch", "Intersection of all species so far");
    for (arg = 1; (arg < argc) &&
		  (count = bl_gff_region_load(&r1, argv[arg])) == 0; ++arg)
	printf("%-20s %2c %2c %2c\n", BL_GFF_REGION_SPECIES(&r1),
		'-', '-', '-');
    if ( arg == argc )
	return EX_OK;
    printf("%-20s %2zu %2c %2s ", BL_GFF_REGION_SPECIES(&r1),
	   BL_GFF_REGION_COUNT(&r1), '*',
	   BL_GFF_REGION_CHROM(&r1));
    print_region_feature_names(&r1);
    putchar('\n');
    max_count = count;

    for (++arg; (arg < argc) && (strcmp(argv[arg], "--diverged") != 0) &&
		  (count = bl_gff_region_load(&rn, argv[arg])) == 0;
		  ++arg)
	printf("%-20s %2c %2c %2c\n", BL_GFF_REGION_SPECIES(&rn),
		'-', '-', '-');
    if ( count > max_count )
	max_count = count;
    if ( arg == argc )
	return EX_OK;
    
    intersect = bl_gff_region_intersect(&r1, &rn);
    printf("%-20s %2zu %2zu %2s ", BL_GFF_REGION_SPECIES(&rn),
	   BL_GFF_REGION_COUNT(&rn), BL_GFF_REGION_COUNT(intersect),
	   BL_GFF_REGION_CHROM(&rn));
    print_region_feature_names(intersect);
    putchar('\n');
    
    snprintf(intersect_file, PATH_MAX + 1, "Intersects/%s-%s-%s",
	     BL_GFF_REGION_GOI(&r1), BL_GFF_REGION_SPECIES(&r1),
	     BL_GFF_REGION_SPECIES(&rn));
    new_intersect_count = BL_GFF_REGION_COUNT(intersect);
    
    for (++arg; (arg < argc) && (strcmp(argv[arg], "--diverged") != 0); ++arg)
    {
	if ( (count = bl_gff_region_load(&rn, argv[arg])) == 0 )
	    printf("%-20s %2c %2c %2c\n", BL_GFF_REGION_SPECIES(&rn),
		    '-', '-', '-');
	else
	{
	    old_intersect_count = new_intersect_count; //BL_GFF_REGION_COUNT(intersect);
	    new_intersect = bl_gff_region_intersect(intersect, &rn);
	    new_intersect_count = BL_GFF_REGION_COUNT(new_intersect);
	    //fprintf(stderr, "new_intersect_count = %d\n", new_intersect_count);
	    if ( count > max_count )
		max_count = count;
	    bl_gff_region_free(intersect);
	    free(intersect);
	    intersect = new_intersect;
	    
	    printf("%-20s %2zu %2zu %2s ", BL_GFF_REGION_SPECIES(&rn),
		    BL_GFF_REGION_COUNT(&rn),
		    BL_GFF_REGION_COUNT(intersect),
		    BL_GFF_REGION_CHROM(&rn));
	    print_region_feature_names(intersect);
	    putchar('\n');
	    
	    if ( new_intersect_count > old_intersect_count )
	    {
		fprintf(stderr, "Ooops, new intersect count cannot exceed old!\n");
		exit(EX_SOFTWARE);
	    }

	    // Update GFF filename
	    strlcat(intersect_file, "-", PATH_MAX + 1);
	    strlcat(intersect_file, BL_GFF_REGION_SPECIES(&rn), PATH_MAX + 1);
	    strlcat(intersect_file, ".gff3", PATH_MAX + 1);
	}
    }
    printf("\nGenes: %d  Conserved: %d  Changed: %d\n",
	    max_count, new_intersect_count, max_count - new_intersect_count);
    
    // Write GFF for intersect
    xt_rmkdir("Intersects", 0777);
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

    if ( arg < argc )   // Implies --diverged
    {
	max_count = 0;
	new_intersect_count = 0;
	++arg;
	puts("\nGenes conserved in each species since diverging from the group above:\n");
	for (; arg < argc; ++arg)
	{
	    if ( (count = bl_gff_region_load(&rn, argv[arg])) == 0 )
		printf("%-20s %2c %2c %2c\n", BL_GFF_REGION_SPECIES(&rn),
			'-', '-', '-');
	    else
	    {
		old_intersect_count = BL_GFF_REGION_COUNT(intersect);
		div_intersect = bl_gff_region_intersect(intersect, &rn);
		new_intersect_count = BL_GFF_REGION_COUNT(div_intersect);
		if ( count > max_count )
		    max_count = count;
		
		printf("%-20s %2zu %2zu %2s ", BL_GFF_REGION_SPECIES(&rn),
			BL_GFF_REGION_COUNT(&rn),
			BL_GFF_REGION_COUNT(div_intersect),
			BL_GFF_REGION_CHROM(&rn));
		print_region_feature_names(div_intersect);
		putchar('\n');
		
		if ( new_intersect_count > old_intersect_count )
		{
		    fprintf(stderr, "Ooops, new intersect count cannot exceed old!\n");
		    exit(EX_SOFTWARE);
		}
	    }
	}
	printf("\nGenes: %d  Conserved: %d  Changed: %d\n",
		max_count, new_intersect_count, max_count - new_intersect_count);
    }
    
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
    fprintf(stderr, "Usage: %s region.gff3 region.gff3 [...]\n"
		    "       [--diverged region.gff3 [region.gff3 ...]]\n",
	    argv[0]);
    exit(EX_USAGE);
}
