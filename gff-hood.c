#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>      // PATH_MAX
#include <biolibc/gff.h>
#include <xtend/string.h>   // Linux strlcpy()
#include <xtend/file.h>
#include <xtend/mem.h>
#include "gff-hood.h"

/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff-hood.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Load a gene neighborhood from a GFF3 file, typically generated
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

int     bl_gff_region_load(bl_gff_region_t *hood, const char *filename)

{
    FILE        *infile;
    int         c;
    char        temp_filename[PATH_MAX + 1], *p, *species, *gene_name;
    const char  *basename;
    bl_gff_t    temp_feature;
    
    if ( (infile = xt_fopen(filename, "r")) == NULL )
	return 0;

    // Parse species and gene name from filename
    if ( (basename = strrchr(filename, '/')) == NULL )
	basename = filename;
    else
	++basename;
    strlcpy(temp_filename, basename, PATH_MAX + 1);
    p = temp_filename;
    species = strsep(&p, "-");
    gene_name = strsep(&p, ".");
    //printf("Loading %s %s\n", species, gene_name);
    
    // Free memory in previously used objects
    if ( hood->species != NULL )
    {
	free(hood->species);
	hood->species = NULL;
    }
    if ( hood->goi != NULL )
    {
	free(hood->goi);
	hood->goi = NULL;
    }
    if ( ((hood->species = strdup(species)) == NULL) ||
	 ((hood->goi = strdup(gene_name)) == NULL) )
    {
	fprintf(stderr, "bl_gff_region_load(): Could not allocate strings.\n");
	return 0;
    }
    
    if ( hood->count == 0 )
    {
	hood->count = 16;
	if ( (hood->features = xt_malloc(hood->count,
	      sizeof(*hood->features))) == NULL )
	{
	    fprintf(stderr, "%s: Could not allocate features array.\n", __FUNCTION__);
	    return 0;
	}
    }
    
    bl_gff_skip_header(infile);
    bl_gff_init(&temp_feature);
    for (c = 0; (bl_gff_read(&temp_feature, infile, BL_GFF_FIELD_ALL)
		 == BL_READ_OK); )
    {
	if ( c == hood->count )
	{
	    hood->count *= 2;
	    if ( (hood->features = xt_realloc(hood->features, hood->count,
		  sizeof(*hood->features))) == NULL )
	    {
		fprintf(stderr, "%s: Could not expand hood->features.\n",
			__FUNCTION__);
		return 0;
	    }
	}
	bl_gff_copy(&hood->features[c], &temp_feature);
	/*printf("%s %" PRIu64 "\n",
		BL_GFF_SEQID(&hood->features[c]),
		BL_GFF_START(&hood->features[c]));*/
	if ( strcasecmp(gene_name, BL_GFF_FEATURE_NAME(&temp_feature)) == 0 )
	    hood->goi_index = c;
	++c;
    }
    hood->count = c;
    fclose(infile);
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
 *  2022-02-20  Jason Bacon Begin
 ***************************************************************************/

void    bl_gff_region_init(bl_gff_region_t *hood)

{
    hood->count = 0;
    hood->goi_index = 0;
    hood->features = NULL;
    hood->species = NULL;
    hood->goi = NULL;
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

void    bl_gff_region_destroy(bl_gff_region_t *hood)

{
    if ( hood->features != NULL )
	free(hood->features);
    if ( hood->species != NULL )
	free(hood->species);
    if ( hood->goi != NULL )
	free(hood->goi);
}


/***************************************************************************
 *  Description:
 *      Generate a commonality score for two gene neighborhoods r1 and r2.
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
    return n - 1;   // Don't count GOI
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
    size_t          n = 0, c1, c2;
    char            *n1, *n2;
    bl_gff_region_t *intersect;

    if ( (intersect = xt_malloc(1, sizeof(*intersect))) == NULL )
    {
	fprintf(stderr, "%s: Could not allocate intersect.\n", __FUNCTION__);
	return NULL;
    }
    bl_gff_region_init(intersect);
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
	    if ( (strcasecmp(n1, n2) == 0) && (strcmp(n1, "unnamed") != 0) &&
		 (strcasecmp(n1, intersect->goi) != 0) )
	    {
		if ( n == intersect->count )
		{
		    if ( (intersect->features = xt_realloc(intersect->features,
			  1, sizeof(*intersect->features))) == NULL )
		    {
			fprintf(stderr, "%s: Could not expand features.\n", __FUNCTION__);
			free(intersect->goi);
			free(intersect->species);
			free(intersect);
			return NULL;
		    }
		}
		
		/*
		 *  Only the feature ID and feature name are common to
		 *  both species and only the name is useful for comparison,
		 *  so leave other fields blank
		 */
		bl_gff_init(&intersect->features[n]);
		bl_gff_set_feature_name(&intersect->features[n], strdup(n2));
		++n;
	    }
	}
    }
    intersect->count = n;
    //fprintf(stderr, "Returning %zu\n", n);
    return intersect;
}
