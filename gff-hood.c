#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>      // PATH_MAX
#include <biolibc/gff.h>
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
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-16  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_hood_load(bl_gff_hood_t *hood, const char *filename)

{
    FILE        *infile;
    int         c;
    char        temp_filename[PATH_MAX + 1], *p, *species, *gene_name;
    const char  *basename;
    
    if ( (infile = fopen(filename, "r")) == NULL )
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
    
    if ( ((hood->species = strdup(species)) == NULL) ||
	 ((hood->goi = strdup(gene_name)) == NULL) )
    {
	fprintf(stderr, "bl_gff_hood_load(): Could not allocate strings.\n");
	return 0;
    }
    
    // FIXME: Maybe hash IDs or names ahead of time to eliminate strcmp()
    // in other functions
    bl_gff_skip_header(infile);
    for (c = 0; (c < BL_GFF_HOOD_MAX_FEATURES) &&
	    (bl_gff_read(&hood->features[c], infile,
		BL_GFF_FIELD_ALL) == BL_READ_OK); )
    {
	/*printf("%s %" PRIu64 "\n",
		BL_GFF_SEQID(&hood->features[c]),
		BL_GFF_START(&hood->features[c]));*/
	if ( strcasecmp(gene_name, BL_GFF_FEATURE_NAME(&hood->features[c])) == 0 )
	    hood->goi_index = c;
	++c;
    }
    hood->gene_count = c;
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

void    bl_gff_hood_init(bl_gff_hood_t *hood)

{
    hood->gene_count = 0;
    hood->goi_index = 0;
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

void    bl_gff_hood_free(bl_gff_hood_t *hood)

{
    if ( hood->species != NULL )
	free(hood->species);
    if ( hood->goi != NULL )
	free(hood->goi);
}


/***************************************************************************
 *  Description:
 *      Generate a commonality score for two gene neighborhoods h1 and h2.
 *      The score represents how many genes they have in common and in
 *      the same order.
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-16  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_hood_commonality(bl_gff_hood_t *h1, bl_gff_hood_t *h2)

{
    int     n = 0, c1, c2;
    char    *n1, *n2;
    
    for (c1 = 0; c1 < h1->gene_count; ++c1)
    {
	n1 = BL_GFF_FEATURE_NAME(&h1->features[c1]);
	for (c2 = 0; c2 < h2->gene_count; ++c2)
	{
	    n2 = BL_GFF_FEATURE_NAME(&h2->features[c2]);
	    if ( (strcasecmp(n1, n2) == 0) && (strcmp(n1, "unnamed") != 0) )
		++n;
	}
    }
    return n - 1;   // Don't count GOI
}
