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

#define MAX_FEATURES    64

int     load_neighborhood(bl_gff_t features[], const char *filename);
int     neighborhood_commonality(bl_gff_t *h1, bl_gff_t *h2, int feature_count);
void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    int         arg, old, new, r, c, old_count, new_count;
    // C is row-major, so features[n] is a contiguous array of MAX_FEATURES
    bl_gff_t    features[2][MAX_FEATURES];
    
    if ( argc < 3 )
	usage(argv);
    
    for (r = 0; r <= 1; ++r)
	for (c = 0; c < MAX_FEATURES; ++c)
	    bl_gff_init(&features[r][c]);
    
    old = 0, new = 1;
    if ( (old_count = load_neighborhood(features[old], argv[1])) == 0 )
    {
	fprintf(stderr, "%s: Error reading %s: %s\n", argv[0],
		argv[1], strerror(errno));
	return EX_NOINPUT;
    }
    // Toggle old and new between 0 and 1, 1 and 0
    for (arg = 2; arg < argc; ++arg, old = (old + 1) % 2, new = (new + 1) % 2)
    {
	printf("%d %d\n", old, new);
	if ( (new_count = load_neighborhood(features[new], argv[arg])) == 0 )
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
	
	printf("%s %s %d\n", argv[arg-1], argv[arg],
	    neighborhood_commonality(features[old], features[new], old_count));
    }
    
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Load a gene neighborhood from a GFF3 file, typically generated
 *      by ms-extract.  The file should contain a short list of adjacent
 *      gene features and nothing more.
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-16  Jason Bacon Begin
 ***************************************************************************/

int     load_neighborhood(bl_gff_t features[], const char *filename)

{
    FILE        *infile;
    int         c;
    
    if ( (infile = fopen(filename, "r")) == NULL )
	return 0;
    
    bl_gff_skip_header(infile);

    c = 0;
    while ( bl_gff_read(&features[c], infile, BL_GFF_FIELD_ALL) == BL_READ_OK )
    {
	printf("%s %" PRIu64 "\n",
		BL_GFF_SEQID(&features[c]), BL_GFF_START(&features[c]));
	++c;
    }
    fclose(infile);
    return c;
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

int     neighborhood_commonality(bl_gff_t *h1, bl_gff_t *h2, int count)

{
    int     n = 0, c;
    char    *n1, *n2;
    
    for (c = 0; c < count; ++c)
    {
	// Just a stub.  This needs to be much smarter.
	n1 = BL_GFF_FEATURE_NAME(&h1[c]);
	n2 = BL_GFF_FEATURE_NAME(&h2[c]);
	if ( (n1 != NULL) && (n2 != NULL) && (strcmp(n1, n2) == 0) )
	    ++n;
    }
    return n;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s neighborhood1.gff3 neighborhood2.gff3 [...]\n",
	    argv[0]);
    exit(EX_USAGE);
}
