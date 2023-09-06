/***************************************************************************
 *  Description:
 *      Add orthologous gene names to a gene list.
 *
 *      The gene list file is a plain text file with one gene per line.
 *
 *      The orthologs file is a TSV file with a gene name from the
 *      species of interest in column 1 and orthologous gene names
 *      from other species in subsequent columns.  Columns may be blank.
 *      This file may be generated by
 *      https://www.ensembl.org/biomart/martview/.
 *      See Genes/add-ortho.sh for more details.
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

#define GENE_ARRAY_SIZE     128
#define MAX_PRINTED_STRINGS 1024

void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    char    *gene_file, *ortho_file,
	    gene[GENE_ARRAY_SIZE], *ortho_gene, *alt_gene,
	    *printed_strings[MAX_PRINTED_STRINGS];
    FILE    *gene_stream, *ortho_stream;
    xt_dsv_line_t  *ortho_line = xt_dsv_line_new();
    int     cmp_status, delim, pc, c, c2;
    
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
     *  leap frog through the two files, matching primary species gene names
     *  (the only gene names in gene_file and the first column in ortho_file).
     */
    
    // Skip header in ortho file if present and read first data line
    xt_dsv_line_init(ortho_line);
    if ( (delim = xt_dsv_line_read(ortho_line, ortho_stream, "\t")) == EOF )
    {
	fprintf(stderr, "add-ortho: Unexpected EOF reading first line of orthologs.\n");
	return EX_DATAERR;
    }
    if ( strcmp(xt_dsv_line_get_fields_ae(ortho_line, 0), "Gene name") == 0 )
    {
	if ( (delim = xt_dsv_line_read(ortho_line, ortho_stream, "\t")) == EOF )
	{
	    fprintf(stderr, "add-ortho: Unexpected EOF reading second line of orthologs.\n");
	    return EX_DATAERR;
	}
    }
    ortho_gene = xt_dsv_line_get_fields_ae(ortho_line, 0);

    while ( xt_fgetline(gene_stream, gene, GENE_ARRAY_SIZE) != EOF )
    {
	pc = 0;
	printf("%s", gene);
	if ( (printed_strings[pc++] = strdup(gene)) == NULL )
	{
	    fprintf(stderr, "Could not allocate printed_strings[pc].\n");
	    return EX_UNAVAILABLE;
	}

	// Skip genes in ortholog list lexically less than GOI
	cmp_status = strcasecmp(ortho_gene, gene);
	while ( (cmp_status < 0) && (delim != EOF ) )
	{
	    xt_dsv_line_free(ortho_line);
	    delim = xt_dsv_line_read(ortho_line, ortho_stream, "\t");
	    ortho_gene = xt_dsv_line_get_fields_ae(ortho_line, 0);
	    cmp_status = strcasecmp(ortho_gene, gene);
	}
	
	while ( (cmp_status == 0) && (delim != EOF) )
	{
	    // Print all uniquely named orthologs
	    for (c = 1; c < xt_dsv_line_get_num_fields(ortho_line); ++c)
	    {
		alt_gene = xt_dsv_line_get_fields_ae(ortho_line, c);
		if ( !xt_strblank(alt_gene) )
		{
		    bool dup = false;
		    // Compare to all other alt genes so far
		    for (c2 = 0; (c2 < pc) && ! dup; ++c2)
			if ( strcasecmp(printed_strings[c2], alt_gene) == 0 )
			    dup = true;
		    if ( ! dup )
		    {
			printf("|%s", alt_gene);
			printed_strings[pc++] = strdup(alt_gene);
		    }
		}
	    }
	    
	    // In case of duplicate GOI lines in the ortholog file
	    xt_dsv_line_free(ortho_line);
	    delim = xt_dsv_line_read(ortho_line, ortho_stream, "\t");
	    ortho_gene = xt_dsv_line_get_fields_ae(ortho_line, 0);
	    cmp_status = strcasecmp(ortho_gene, gene);
	}
	putchar('\n');
	for (c2 = 0; c2 < pc; ++c2)
	    free(printed_strings[c2]);
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
