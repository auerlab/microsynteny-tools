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
#include <xtend/mem.h>
#include <xtend/file.h>
#include "ms-extract.h"

int     main(int argc,char *argv[])

{
    bl_gff_t    feature = BL_GFF_INIT;
    FILE        *gff_stream, *header_stream, *gene_stream;
    char        *gene_file,
		**gene_names,
		*gff_file,
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
	    adjacent_genes = strtoul(argv[++arg], &end, 10);
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
    
    if ( arg != argc - 2 )
	usage(argv);
    gff_file = argv[arg++];
    gene_file = argv[arg];

    /*
     *  GFF file must be seekable so we can back up and extrac the upstream
     *  genes, so it must be uncompressed.  Don't use xt_fopen() here, which
     *  will call popen() to pipe compressed data in.
     */
    if ( (gff_stream = fopen(gff_file, "r")) == NULL )
    {
	fprintf(stderr, "ms-extract: Could not open %s: %s\n",
		gff_file, strerror(errno));
	return EX_NOINPUT;
    }
    header_stream = bl_gff_skip_header(gff_stream);
    
    if ( (gene_stream = xt_fopen(gene_file, "r")) == NULL )
    {
	fprintf(stderr, "ms-extract: Could not open %s: %s\n",
		gene_file, strerror(errno));
	return EX_NOINPUT;
    }
    if ( (gene_count = xt_inhale_strings(gene_stream, &gene_names)) == 0 )
    {
	fprintf(stderr, "ms-extract: Unable to read genes from %s.\n", gene_file);
	return EX_NOINPUT;
    }
    fclose(gene_stream);
    
    if ( (gff_basename = strrchr(gff_file, '/')) == NULL )
	gff_basename = gff_file;
    else
	++gff_basename; // First char after /
    
    // Note: Destroys gff_file == argv[1]
    if ( (ext = strchr(gff_basename, '.')) != NULL )
	*ext = '\0';

    // For consistent output and easy comparison by other programs
    for (c = 0; c < gene_count; ++c)
    {
	//fprintf(stderr, "gene: %s\n", gene_names[c]);
	strlower(gene_names[c]);
    }
    
    genes_found = 0;
    while ( (genes_found < gene_count) &&
	    (bl_gff_read(&feature, gff_stream, BL_GFF_FIELD_ALL) == BL_READ_OK) )
    {
	// Only interested in gene features
	if ( strcasecmp(BL_GFF_TYPE(&feature), "gene") == 0 )
	{
	    /*
	     *  Index positions of all genes in the file so we can
	     *  quickly seek() them later.
	     */
	    if ( bl_gff_index_add(&gi, &feature) != BL_GFF_INDEX_OK )
	    {
		fprintf(stderr, "bl_gff_index_add_pos() failed.\n");
		return EX_SOFTWARE;
	    }
	    //printf("Added %s to index.\n", BL_GFF_FEATURE_NAME(&feature));
	    
	    for (c = 0; c < gene_count; ++c)
	    {
		// Ensembl GFFs are not consistent with capitalization.
		// All lower case for some species, capitalized for others
		if ( strcasecmp(BL_GFF_FEATURE_NAME(&feature), gene_names[c]) == 0 )
		{
		    ++genes_found;
		    printf("\n%s %s:\n", gff_basename, gene_names[c]);
		    
		    // Path name format is important, parsed by ms-divergence
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
		    break;  // Don't bother comparing the other genes
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

int     extract_neighborhood(bl_gff_t *goi, bl_gff_index_t *gi,
	    FILE *gff_stream, FILE *header_stream, char *hood_file,
	    uint64_t adjacent_genes, uint64_t max_nt_distance)

{
    bl_gff_t    neighbor = BL_GFF_INIT;
    char        *neighbor_name;
    uint64_t    g, len;
    int64_t     distance;
    FILE        *hood_stream;
    
    if ( (hood_stream = fopen(hood_file, "w")) == NULL )
    {
	fprintf(stderr, "ms-extract: Cannot open %s: %s\n",
		hood_file, strerror(errno));
	return EX_CANTCREAT;
    }
    // Add the original GFF header to each neighborhood GFF
    bl_gff_copy_header(header_stream, hood_stream);
    
    /*
     *  Back up to at the position of this gene - adjacent_genes, but no
     *  more than max_nt_distance nucleotides, and output all genes from
     *  there to position + (adjacent_genes or max_nt_distance)
     */
    if ( bl_gff_index_seek_reverse(gi, gff_stream, goi,
	    adjacent_genes, max_nt_distance) != 0 )
    {
	fprintf(stderr, "ms-extract: Seek %zu failed.\n",
		BL_GFF_FILE_POS(goi));
	return EX_SOFTWARE;
    }
    
    printf("%-3s %10s %10s %10s %10s %s\n",
	    "Chr", "Start", "End", "Len", "Distance", "Name");
    // From leftmost neighbor read adjacent_genes before and after GOI
    for (g = 0; (g < adjacent_genes * 2 + 1) &&
		bl_gff_read(&neighbor, gff_stream,
			BL_GFF_FIELD_ALL) == BL_READ_OK; )
    {
	if ( strcmp(BL_GFF_TYPE(&neighbor), "gene") == 0 )
	{
	    // Genes after GOI have not been indexed yet
	    if ( g > adjacent_genes )
	    {
		if ( bl_gff_index_add(gi, &neighbor) != BL_GFF_INDEX_OK )
		{
		    fprintf(stderr, "bl_gff_index_add_pos() failed.\n");
		    return EX_SOFTWARE;
		}
		//printf("Added %s to index.\n", BL_GFF_FEATURE_NAME(&neighbor));
	    }
	    ++g;
	    neighbor_name = BL_GFF_FEATURE_NAME(&neighbor);
	    if ( neighbor_name == NULL )
		neighbor_name = "unnamed";
	    if ( BL_GFF_START(&neighbor) < BL_GFF_START(goi) )
		distance = BL_GFF_END(&neighbor) - BL_GFF_START(goi);
	    else if ( BL_GFF_START(&neighbor) > BL_GFF_START(goi) )
		distance = BL_GFF_START(&neighbor) - BL_GFF_END(goi);
	    else
		distance = 0;
	    len = BL_GFF_END(&neighbor) - BL_GFF_START(&neighbor);
	    
	    printf("%-3s %10" PRIu64 " %10" PRIu64 " %10" PRId64 " %10" PRId64 " %s\n",
		    BL_GFF_SEQID(&neighbor),
		    BL_GFF_START(&neighbor),
		    BL_GFF_END(&neighbor), len, distance, neighbor_name);
	
	    // Hack gene name into source field for DNA Feature Viewer.
	    // There should be a way to avoid this.
	    bl_gff_set_source_cpy(&neighbor, neighbor_name,
		BL_GFF_SOURCE_MAX_CHARS);
	    
	    // Output gene neighborhood in GFF format for above
	    bl_gff_write(&neighbor, hood_stream,
		BL_GFF_FIELD_ALL);
	    
	    // FIXME: Restructure the loop so this is a loop condition
	    if ( BL_GFF_START(&neighbor) >
		    BL_GFF_END(goi) + max_nt_distance )
		break;
	}
	bl_gff_free(&neighbor);
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
		    "   species.gff3 gene-list.txt\n", argv[0]);
    exit(EX_USAGE);
}

