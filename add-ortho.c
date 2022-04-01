/***************************************************************************
 *  Description:
 *  
 *  Arguments:
 *
 *  Returns:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-01  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>
#include <xtend/string.h>
#include <xtend/file.h>
#include <xtend/dsv.h>

#define GENE_ARRAY_SIZE 128

void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    char    *gene_file, *ortho_file,
	    gene[GENE_ARRAY_SIZE], *ortho_gene, *alt_gene;
    FILE    *gene_stream, *ortho_stream;
    dsv_line_t  ortho_line;
    int     cmp_status, delim;
    
    switch(argc)
    {
	case 3:
	    gene_file = argv[1];
	    ortho_file = argv[2];
	    break;
	
	default:
	    usage(argv);
    }
    
    if ( (gene_stream = xt_fopen(gene_file, "r")) == NULL )
    {
	fprintf(stderr, "add-ortho: Cannot open %s: %s\n",
		gene_file, strerror(errno));
	return EX_NOINPUT;
    }
    
    if ( (ortho_stream = xt_fopen(ortho_file, "r")) == NULL )
    {
	fprintf(stderr, "add-ortho: Cannot open %s: %s\n",
		ortho_file, strerror(errno));
	return EX_NOINPUT;
    }
    
    /*
     *  Both files must be sorted by gene name.  The loops below will
     *  leap frog through the two files, matching zebrafish gene names
     *  (the only gene names in gene_file and the first column in ortho_file).
     */
    
    // Skip header in ortho file if present
    // dsv_line_init(&ortho_line);
    if ( (delim = dsv_line_read(&ortho_line, ortho_stream, "\t")) == EOF )
    {
	fprintf(stderr, "add-ortho: Unexpected EOF reading first line of orthologs.\n");
	return EX_DATAERR;
    }
    if ( strcmp(DSV_LINE_FIELDS_AE(&ortho_line, 0), "Gene name") == 0 )
    {
	if ( (delim = dsv_line_read(&ortho_line, ortho_stream, "\t")) == EOF )
	{
	    fprintf(stderr, "add-ortho: Unexpected EOF reading second line of orthologs.\n");
	    return EX_DATAERR;
	}
    }
    ortho_gene = DSV_LINE_FIELDS_AE(&ortho_line, 0);
    while ( xt_fgetline(gene_stream, gene, GENE_ARRAY_SIZE) != EOF )
    {
	printf("%s", gene);
	cmp_status = strcasecmp(ortho_gene, gene);
	while ( (cmp_status < 0) && (delim != EOF ) )
	{
	    //printf("\nSkipping %s...", ortho_gene);
	    dsv_line_free(&ortho_line);
	    delim = dsv_line_read(&ortho_line, ortho_stream, "\t");
	    ortho_gene = DSV_LINE_FIELDS_AE(&ortho_line, 0);
	    cmp_status = strcasecmp(ortho_gene, gene);
	}
	//printf("\nStopped at %s.\n", ortho_gene);
	if ( cmp_status == 0 )
	{
	    for (int c = 1; c < DSV_LINE_NUM_FIELDS(&ortho_line); ++c)
	    {
		alt_gene = DSV_LINE_FIELDS_AE(&ortho_line, c);
		if ( !strblank(alt_gene) )
		{
		    bool dup = false;
		    // Compare to all other alt genes so far
		    for (int c2 = 0; (c2 < c) && ! dup; ++c2)
			if ( strcasecmp(DSV_LINE_FIELDS_AE(&ortho_line, c2), alt_gene) == 0 )
			    dup = true;
		    if ( ! dup )
			printf("|%s", alt_gene);
		}
	    }
	}
	putchar('\n');
    }
    
    fclose(ortho_stream);
    fclose(gene_stream);
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s gene-file.txt ortho-file.tsv[.gz|.bz2|.xz]\n",
	    argv[0]);
    exit(EX_USAGE);
}
