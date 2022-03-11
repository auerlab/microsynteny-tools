#ifdef __linux__
#define _GNU_SOURCE         // asprintf()
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>      // PATH_MAX
#include <biolibc/gff.h>
#include <xtend/string.h>   // Linux strlcpy()
#include <xtend/file.h>
#include <xtend/mem.h>
#include "gff-region.h"

/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff-region.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Load a gene region from a GFF3 file, typically generated
 *      by ms-extract.  The file should contain a short list of adjacent
 *      gene features and nothing more.
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
 *  2022-02-16  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_region_load(bl_gff_region_t *region, const char *filename)

{
    FILE        *infile;
    char        temp_filename[PATH_MAX + 1], *p, *species, *gene_name;
    const char  *basename;
    bl_gff_t    temp_feature;

    /*
     *  Set species and goi even if input file does not exist
     */
    
    // Parse species and gene name from filename
    if ( (basename = strrchr(filename, '/')) == NULL )
	basename = filename;
    else
	++basename;
    
    // FIXME: Make a libxtend function for this?
    strlcpy(temp_filename, basename, PATH_MAX + 1);
    p = temp_filename;
    species = strsep(&p, "-");
    gene_name = strsep(&p, "-");
    //printf("Loading %s %s\n", species, gene_name);
    
    // Free memory in previously used objects
    if ( region->species != NULL )
    {
	free(region->species);
	region->species = NULL;
    }
    if ( region->goi != NULL )
    {
	free(region->goi);
	region->goi = NULL;
    }
    if ( region->chrom != NULL )
    {
	free(region->chrom);
	region->chrom = NULL;
    }
    if ( ((region->species = strdup(species)) == NULL) ||
	 ((region->goi = strdup(gene_name)) == NULL) )
    {
	fprintf(stderr, "bl_gff_region_load(): Could not allocate strings.\n");
	return 0;
    }
    
    if ( (infile = xt_fopen(filename, "r")) == NULL )
	return 0;
    
    if ( region->array_size == 0 )
    {
	region->array_size = 16;
	if ( (region->features = xt_malloc(region->array_size,
	      sizeof(*region->features))) == NULL )
	{
	    fprintf(stderr, "%s: Could not allocate features array.\n", __FUNCTION__);
	    return 0;
	}
    }
    
    bl_gff_skip_header(infile);
    bl_gff_init(&temp_feature);
    for (region->count = 0; (bl_gff_read(&temp_feature, infile, BL_GFF_FIELD_ALL)
		 == BL_READ_OK); )
    {
	if ( region->count == region->array_size )
	{
	    region->array_size *= 2;
	    if ( (region->features = xt_realloc(region->features, region->array_size,
		  sizeof(*region->features))) == NULL )
	    {
		fprintf(stderr, "%s: Could not expand region->features.\n",
			__FUNCTION__);
		return 0;
	    }
	}
	bl_gff_copy(&region->features[region->count], &temp_feature);
	/*printf("%s %" PRId64 "\n",
		BL_GFF_SEQID(&region->features[c]),
		BL_GFF_START(&region->features[c]));*/
	//printf("%s %s\n", gene_name, BL_GFF_FEATURE_NAME(&temp_feature));
	if ( strcasecmp(gene_name, BL_GFF_FEATURE_NAME(&temp_feature)) == 0 )
	{
	    region->goi_index = region->count;
	    if ( (region->chrom = strdup(BL_GFF_SEQID(&temp_feature))) == NULL )
	    {
		fprintf(stderr, "bl_gff_region_load(): Failed to strdup() chrom.\n");
		return 0;
	    }
	}
	++region->count;
    }
    fclose(infile);
    return region->count;
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
 *  2022-02-20  Jason Bacon Begin
 ***************************************************************************/

void    bl_gff_region_init(bl_gff_region_t *region)

{
    region->array_size = 0;
    region->count = 0;
    region->goi_index = 0;
    region->features = NULL;
    region->species = NULL;
    region->goi = NULL;
    region->chrom = NULL;
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
 *  2022-02-17  Jason Bacon Begin
 ***************************************************************************/

void    bl_gff_region_free(bl_gff_region_t *region)

{
    if ( region->features != NULL )
	free(region->features);
    if ( region->species != NULL )
	free(region->species);
    if ( region->goi != NULL )
	free(region->goi);
    if ( region->chrom != NULL )
	free(region->chrom);
    bl_gff_region_init(region);
}


/***************************************************************************
 *  Description:
 *      Generate a commonality score for two gene regions r1 and r2.
 *      The score represents how many genes they have in common and in
 *      the same order.
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-16  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_region_commonality(bl_gff_region_t *r1, bl_gff_region_t *r2)

{
    int     n = 0, c1, c2;
    char    *n1, *n2;
    
    for (c1 = 0; c1 < r1->count; ++c1)
    {
	n1 = BL_GFF_FEATURE_NAME(&r1->features[c1]);
	for (c2 = 0; c2 < r2->count; ++c2)
	{
	    n2 = BL_GFF_FEATURE_NAME(&r2->features[c2]);
	    if ( (strcasecmp(n1, n2) == 0) && (strcmp(n1, "unnamed") != 0) )
		++n;
	}
    }
    return n;
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

bl_gff_region_t   *bl_gff_region_intersect(bl_gff_region_t *r1, bl_gff_region_t *r2)

{
    size_t          c1, c2;
    char            *n1, *n2, *attr;
    bl_gff_region_t *intersect;

    if ( (intersect = xt_malloc(1, sizeof(*intersect))) == NULL )
    {
	fprintf(stderr, "%s: Could not allocate intersect.\n", __FUNCTION__);
	return NULL;
    }
    bl_gff_region_init(intersect);
    intersect->array_size = 16;
    if ( (intersect->features = xt_malloc(intersect->array_size,
	  sizeof(*intersect->features))) == NULL )
    {
	fprintf(stderr, "%s: Could not allocate features array.\n", __FUNCTION__);
	return 0;
    }
    
    if ( (intersect->species = strdup("Intersection")) == NULL )
    {
	fprintf(stderr, "%s: Could not strdup species.\n", __FUNCTION__);
	free(intersect);
	return NULL;
    }
    
    if ( (intersect->goi = strdup(r1->goi)) == NULL )
    {
	fprintf(stderr, "%s: Could not strdup GOI.\n", __FUNCTION__);
	free(intersect->species);
	free(intersect);
	return NULL;
    }
    
    //fprintf(stderr, "%s %s\n", r1->species, r2->species);
    // FIXME: Don't count multiple copies of the same gene
    for (c1 = 0; c1 < r1->count; ++c1)
    {
	n1 = BL_GFF_FEATURE_NAME(&r1->features[c1]);
	for (c2 = 0; c2 < r2->count; ++c2)
	{
	    n2 = BL_GFF_FEATURE_NAME(&r2->features[c2]);
	    // fprintf(stderr, "%s %s\n", n1, n2);
	    //     (strcasecmp(n1, intersect->goi) != 0) &&
	    if ( (strcasecmp(n1, n2) == 0) && (strcmp(n1, "unnamed") != 0) &&
		 !bl_gff_region_duplicate_gene(intersect, n1) )
	    {
		if ( intersect->count == intersect->array_size )
		{
		    intersect->array_size *= 2;
		    if ( (intersect->features = xt_realloc(intersect->features,
			  intersect->array_size, sizeof(*intersect->features))) == NULL )
		    {
			fprintf(stderr, "%s: Could not expand features.\n", __FUNCTION__);
			free(intersect->goi);
			free(intersect->species);
			free(intersect);
			return NULL;
		    }
		}
		
		/*
		 *  Only the feature name is common to both species and
		 *  so leave other fields blank.
		 */
		bl_gff_init(&intersect->features[intersect->count]);
		bl_gff_set_type_cpy(&intersect->features[intersect->count],
		    "gene", BL_GFF_TYPE_MAX_CHARS + 1);
		strlower(n2);
		bl_gff_set_feature_name(&intersect->features[intersect->count],
		    strdup(n2));
		asprintf(&attr, "Name=%s;", n2);
		bl_gff_set_attributes(&intersect->features[intersect->count],
		    attr);
		++intersect->count;
	    }
	}
    }
    //fprintf(stderr, "Returning %zu\n", n);
    return intersect;
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
 *  2022-03-02  Jason Bacon Begin
 ***************************************************************************/

bool    bl_gff_region_duplicate_gene(bl_gff_region_t *r, const char *feature_name)

{
    int     c;
    
    /*
     *  FIXME: Maybe find a faster way.  Linear search is fine for now
     *  since we are using very small lists (no more than 20) for
     *  microsynteny-tools.
     */
    
    for (c = 0; c < r->count; ++c)
	if ( strcmp(BL_GFF_FEATURE_NAME(&r->features[c]), feature_name) == 0 )
	{
	    fprintf(stderr, "Duplicate gene: %s\n", feature_name);
	    return true;
	}
    return false;
}
