/***************************************************************************
 *  Description:
 *      Show genes in the neighborhood of the specified gene for
 *      the specified species.
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-01  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/param.h>      // PATH_MAX
#include <biolibc/gff.h>
#include <biolibc/gff-index.h>
#include <xtend/string.h>
#include "ms-extract.h"


int     main(int argc,char *argv[])

{
    bl_gff_t    feature = BL_GFF_INIT;
    FILE        *gff_stream, *header_stream;
    char        **gene_names,
		*gff_filename,
		*end,
		*output_dir = ".",
		*gff_basename,
		*ext,
		hood_file[PATH_MAX];
    uint64_t    max_nt_distance = DEFAULT_NT_DISTANCE,
		adjacent_genes = DEFAULT_ADJACENT_GENES;
    int         arg,
		c,
		status,
		gene_count,
		genes_found;
    bl_gff_index_t    gi = BL_GFF_INDEX_INIT;

    for (arg = 1; (arg < argc) && (argv[arg][0] == '-'); ++arg)
    {
	if ( strcmp(argv[arg], "--output-dir") == 0 )
	    output_dir = argv[++arg];
	else if ( strcmp(argv[arg], "--max-nt-distance") == 0 )
	{
	    max_nt_distance = strtoul(argv[++arg], &end, 10);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "ms-extract: Invalid --max-nt-distance: %s\n",
			argv[arg]);
		usage(argv);
	    }
	}
	else if ( strcmp(argv[arg], "--adjacent-genes") == 0 )
	{
	    gene_count = strtoul(argv[++arg], &end, 10);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "ms-extract: Invalid --adjacent-genes: %s\n",
			argv[arg]);
		usage(argv);
	    }
	}
	else
	    usage(argv);
    }
    
    if ( arg > argc - 2 )
	usage(argv);
    gff_filename = argv[arg++];
    gene_names = argv + arg;
    gene_count = argc - arg;
    
    // FIXME: Maybe make xt_fopen_seekable() function to automatically
    // unzip compressed files if present
    if ( (gff_stream = fopen(gff_filename, "r")) == NULL )
    {
	fprintf(stderr, "ms-extract: Could not open %s: %s\n",
		gff_filename, strerror(errno));
	exit(EX_NOINPUT);
    }
    header_stream = bl_gff_skip_header(gff_stream);
    
    if ( (gff_basename = strrchr(gff_filename, '/')) == NULL )
	gff_basename = gff_filename;
    // Note: Destroys gff_filename == argv[1]
    if ( (ext = strchr(gff_basename, '.')) != NULL )
	*ext = '\0';

    genes_found = 0;
    while ( (genes_found < gene_count) &&
	    (bl_gff_read(&feature, gff_stream, BL_GFF_FIELD_ALL) == BL_READ_OK) )
    {
	if ( strcasecmp(BL_GFF_TYPE(&feature), "gene") == 0 )
	{
	    // Index positions of all genes in the file
	    if ( bl_gff_index_add(&gi, &feature) != BL_GFF_INDEX_OK )
	    {
		fprintf(stderr, "bl_gff_index_add_pos() failed.\n");
		return EX_SOFTWARE;
	    }
	    
	    for (c = 0; c < gene_count; ++c)
	    {
		if ( (strcasecmp(BL_GFF_FEATURE_NAME(&feature), gene_names[c]) == 0) )
		{
		    ++genes_found;
		    printf("\nFound %s:\n", gene_names[c]);
		    strlower(gene_names[c]);
		    // FIXME: Factor out to a bl_gff_hood function
		    // Name format is important, used by ms-divergence
		    snprintf(hood_file, PATH_MAX, "%s/%s-%s.gff3",
			output_dir, gff_basename, gene_names[c]);
		    status = extract_neighborhood(&feature, &gi, gff_stream,
			header_stream, hood_file,
			adjacent_genes, max_nt_distance);
		    if ( status != EX_OK )
		    {
			fclose(gff_stream);
			exit(status);
		    }
		}
	    }
	}
	bl_gff_free(&feature);
    }
    fclose(gff_stream);
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Extract the gene neighborhood of the given feature, aided by
 *      index gi for speed, saving the neighborhood in GFF3 format to
 *      hood_file.
 *
 *  Returns:
 *      EX_OK on success or another sysexits constant on error
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-13  Jason Bacon Begin
 ***************************************************************************/

int     extract_neighborhood(bl_gff_t *feature, bl_gff_index_t *gi,
	    FILE *gff_stream, FILE *header_stream, char *hood_file,
	    uint64_t adjacent_genes, uint64_t max_nt_distance)

{
    bl_gff_t    neighbor_gene = BL_GFF_INIT;
    char        *neighbor_name;
    uint64_t    g;
    FILE        *hood_stream;
    
    if ( (hood_stream = fopen(hood_file, "w")) == NULL )
    {
	fprintf(stderr, "ms-extract: Cannot open %s: %s\n",
		hood_file, strerror(errno));
	return EX_CANTCREAT;
    }
    bl_gff_copy_header(header_stream, hood_stream);
    
    /*
     *  Back up to at the position of this gene - max_nt_distance
     *  and output all genes from there to position + max_nt_distance
     */
    if ( bl_gff_index_seek_reverse(gi, gff_stream, feature,
	    adjacent_genes, max_nt_distance) != 0 )
    {
	fprintf(stderr, "ms-extract: Seek %zu failed.\n",
		BL_GFF_FILE_POS(feature));
	return EX_SOFTWARE;
    }
    
    // From leftmost neighbor read adjacent_genes before and after ref
    for (g = 0; (g < adjacent_genes * 2 + 1) &&
		bl_gff_read(&neighbor_gene, gff_stream,
			BL_GFF_FIELD_ALL) == BL_READ_OK; )
    {
	if ( strcmp(BL_GFF_TYPE(&neighbor_gene), "gene") == 0 )
	{
	    ++g;
	    neighbor_name = BL_GFF_FEATURE_NAME(&neighbor_gene);
	    if ( neighbor_name == NULL )
		neighbor_name = "unnamed";
	    printf("%s\t%" PRIu64 "\t%" PRIu64 "\t%s\n",
		    BL_GFF_SEQID(&neighbor_gene),
		    BL_GFF_START(&neighbor_gene),
		    BL_GFF_END(&neighbor_gene), neighbor_name);
	
	    // Hack gene name into source field for
	    // DNA Feature Viewer.  Should be a way to avoid this.
	    bl_gff_set_source_cpy(&neighbor_gene, neighbor_name,
		BL_GFF_SOURCE_MAX_CHARS);
	    // Output gene neighborhood in GFF format for above
	    bl_gff_write(&neighbor_gene, hood_stream,
		BL_GFF_FIELD_ALL);
	    
	    if ( BL_GFF_START(&neighbor_gene) >
		    BL_GFF_END(feature) + max_nt_distance )
		break;
	}
	bl_gff_free(&neighbor_gene);
    }
    fclose(hood_stream);
    
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s\n"
		    "   [--output-dir dir]\n"
		    "   [--adjacent-genes genes]\n"
		    "   [--max-nt-distance nucleotides]\n"
		    "   file.gff3 gene-name [gene-name ...]\n", argv[0]);
    exit(EX_USAGE);
}

