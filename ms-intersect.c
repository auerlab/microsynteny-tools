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
#include <biolibc/gff3.h>
#include "gff-region.h"

#define SPECIES_MAX_CHARS   256

int     intersect(int argc, char *argv[]);
void    print_region_feature_names(bl_gff3_region_t *region);
void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    int     status;
    
    if ( argc < 4 )
	usage(argv);

    status = intersect(argc, argv);
    
    puts("\n===================================================================");
    puts("Note: This is a conservative estimate of conserved genes.");
    puts("There may be unnamed or analogous genes that are not shown here.");
    puts("Use ms-stack to see a full list of neighbors that might be excluded");
    puts("From the output above.");
    puts("===================================================================");
    
    return status;
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
    int             arg, count, new_intersect_count,
		    max_count;
    bl_gff3_region_t r1, rn, *intersect, *new_intersect, *div_intersect;
    char            // intersect_file[PATH_MAX + 1],
		    previous_species[SPECIES_MAX_CHARS + 1] = "";
    //FILE            *intersect_stream;
    
    if ( argc < 3 )
	usage(argv);
    
    bl_gff3_region_init(&r1);
    bl_gff3_region_init(&rn);
    // output_dir = argv[1];

    // Output header
    printf("%-20s %2s %2s %2s  %s\n", "Species", "Ne", "Co",
	    "Ch", "Intersection of all species so far");
    
    // Skip non-existent region files and load first existent
    for (arg = 2, count = 0; (arg < argc) &&
		  (count = bl_gff3_region_load(&r1, argv[arg])) == 0; ++arg)
	printf("%-20s %2c %2c %2c\n", BL_GFF3_REGION_SPECIES(&r1),
		'-', '-', '-');
    if ( arg == argc )
	return EX_OK;
    
    // Output first existing region
    printf("%-20s %2zu %2c %2s ", BL_GFF3_REGION_SPECIES(&r1),
	   BL_GFF3_REGION_COUNT(&r1), '*',
	   BL_GFF3_REGION_CHROM(&r1));
    print_region_feature_names(&r1);
    putchar('\n');
    max_count = count;

    // Skip non-existent files after first existent
    for (++arg; (arg < argc) && (strcmp(argv[arg], "--diverged") != 0) &&
		  (count = bl_gff3_region_load(&rn, argv[arg])) == 0;
		  ++arg)
	printf("%-20s %2c %2c %2c\n", BL_GFF3_REGION_SPECIES(&rn),
		'-', '-', '-');
    if ( count > max_count )
	max_count = count;
    
    // If no more GFF files, we're done
    if ( arg == argc )
	return EX_OK;
    
    if ( strcmp(argv[arg], "--diverged") == 0 )
    {
	intersect = &r1;
    }
    else
    {
	intersect = bl_gff3_region_intersect(&r1, &rn);
	printf("%-20s %2zu %2zu %2s ", BL_GFF3_REGION_SPECIES(&rn),
	       BL_GFF3_REGION_COUNT(&rn), BL_GFF3_REGION_COUNT(intersect),
	       BL_GFF3_REGION_CHROM(&rn));
	print_region_feature_names(intersect);
	putchar('\n');
	
	/*
	snprintf(intersect_file, PATH_MAX + 1, "%s/%s-%s-%s", output_dir,
		 BL_GFF3_REGION_GOI(&r1), BL_GFF3_REGION_SPECIES(&r1),
		 BL_GFF3_REGION_SPECIES(&rn));
	*/
	new_intersect_count = BL_GFF3_REGION_COUNT(intersect);
	
	for (++arg; (arg < argc) && (strcmp(argv[arg], "--diverged") != 0); ++arg)
	{
	    if ( (count = bl_gff3_region_load(&rn, argv[arg])) == 0 )
		printf("%-20s %2c %2c %2c\n", BL_GFF3_REGION_SPECIES(&rn),
			'-', '-', '-');
	    else
	    {
		//old_intersect_count = new_intersect_count; //BL_GFF3_REGION_COUNT(intersect);
		new_intersect = bl_gff3_region_intersect(intersect, &rn);
		new_intersect_count = BL_GFF3_REGION_COUNT(new_intersect);
		//fprintf(stderr, "new_intersect_count = %d\n", new_intersect_count);
		if ( count > max_count )
		    max_count = count;
		bl_gff3_region_free(intersect);
		free(intersect);
		intersect = new_intersect;
		
		printf("%-20s %2zu %2zu %2s ", BL_GFF3_REGION_SPECIES(&rn),
			BL_GFF3_REGION_COUNT(&rn),
			BL_GFF3_REGION_COUNT(intersect),
			BL_GFF3_REGION_CHROM(&rn));
		print_region_feature_names(intersect);
		putchar('\n');
		
		/* Actually count CAN increase due to GOI orthologs
		if ( new_intersect_count > old_intersect_count )
		{
		    fprintf(stderr, "Ooops1, new intersect count cannot exceed old!\n");
		    exit(EX_SOFTWARE);
		}
		*/
    
		// Update GFF filename
		if ( strcmp(BL_GFF3_REGION_SPECIES(&rn), previous_species) != 0 )
		{
		    //strlcat(intersect_file, "-", PATH_MAX + 1);
		    //strlcat(intersect_file, BL_GFF3_REGION_SPECIES(&rn), PATH_MAX + 1);
		    strlcpy(previous_species, BL_GFF3_REGION_SPECIES(&rn), SPECIES_MAX_CHARS + 1);
		}
	    }
	}
	
	//strlcat(intersect_file, ".gff3", PATH_MAX + 1);
	//fprintf(stderr, "intersect_file = %s\n", intersect_file);
	
	printf("\nGenes: %d  Conserved: %d  Changed: %d\n",
		max_count, new_intersect_count, max_count - new_intersect_count);
	
	// Write GFF for intersect
	/*
	xt_rmkdir(output_dir, 0777);
	if ( (intersect_stream = fopen(intersect_file, "w")) == NULL )
	{
	    fprintf(stderr, "ms-intersect: Could not open %s for write: %s.\n",
		    intersect_file, strerror(errno));
	    exit(EX_CANTCREAT);
	}
	fprintf(intersect_stream, "##gff-version 3\n");
	fprintf(intersect_stream, "## %s\n", intersect_file);
	fprintf(intersect_stream, "## This file is simply a list of genes shared by multiple species.\n");
	fprintf(intersect_stream, "## There is no chrom, start, end, etc. since they will differ across species.\n");
	for (c = 0; c < BL_GFF3_REGION_COUNT(intersect); ++c)
	    bl_gff3_write(&BL_GFF3_REGION_FEATURES_AE(intersect, c),
			 intersect_stream, BL_GFF3_FIELD_ALL);
	fclose(intersect_stream);
	*/
    }
    
    if ( arg < argc )   // Implies --diverged
    {
	max_count = 0;
	new_intersect_count = 0;
	++arg;
	
	puts("\nGenes conserved in each species since diverging from the group above:\n");
	for (; arg < argc; ++arg)
	{
	    if ( (count = bl_gff3_region_load(&rn, argv[arg])) == 0 )
		printf("%-20s %2c %2c %2c\n", BL_GFF3_REGION_SPECIES(&rn),
			'-', '-', '-');
	    else
	    {
		//old_intersect_count = BL_GFF3_REGION_COUNT(intersect);
		div_intersect = bl_gff3_region_intersect(intersect, &rn);
		new_intersect_count = BL_GFF3_REGION_COUNT(div_intersect);
		if ( count > max_count )
		    max_count = count;
		
		printf("%-20s %2zu %2zu %2s ", BL_GFF3_REGION_SPECIES(&rn),
			BL_GFF3_REGION_COUNT(&rn),
			BL_GFF3_REGION_COUNT(div_intersect),
			BL_GFF3_REGION_CHROM(&rn));
		print_region_feature_names(div_intersect);
		putchar('\n');
		
		/* Actually count CAN increase due to GOI orthologs
		if ( new_intersect_count > old_intersect_count )
		{
		    fprintf(stderr, "Ooops2, new intersect count cannot exceed old!\n");
		    exit(EX_SOFTWARE);
		}
		*/
	    }
	}
	//printf("\nGenes: %d  Conserved: %d  Changed: %d\n",
	//        max_count, new_intersect_count, max_count - new_intersect_count);
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

void    print_region_feature_names(bl_gff3_region_t *region)

{
    int         c;
    bl_gff3_t    *gff;
    
    for (c = 0; c < BL_GFF3_REGION_COUNT(region); ++c)
    {
	gff = &(BL_GFF3_REGION_FEATURES_AE(region, c));
	if ( BL_GFF3_STRAND(gff) == '+' )
	    printf(" %s+", BL_GFF3_FEATURE_NAME(gff));
	else if ( BL_GFF3_STRAND(gff) == '-' )
	    printf(" -%s", BL_GFF3_FEATURE_NAME(gff));
	else
	    printf(" %s", BL_GFF3_FEATURE_NAME(gff));
    }
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s output-dir region.gff3 region.gff3 [...]\n"
		    "       [--diverged region.gff3 [region.gff3 ...]]\n",
	    argv[0]);
    exit(EX_USAGE);
}
