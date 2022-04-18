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
#include "alt-str.h"

int     main(int argc,char *argv[])

{
    bl_gff_t    feature;
    FILE        *gff_stream, *header_stream, *gene_stream;
    char        *gene_file, *gff_file,
		*end,
		*output_dir = ".", *gff_basename, *ext,
		region_file[PATH_MAX];
    int64_t     max_nt_distance = DEFAULT_NT_DISTANCE;
    int         adjacent_genes = DEFAULT_ADJACENT_GENES,
		arg, c, status, gene_count;
    bl_gff_index_t    gi = BL_GFF_INDEX_INIT;
    alt_str_t   *gene_names;

    bl_gff_init(&feature);
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
    if ( (gene_count = xt_alt_str_inhale_list(&gene_names, gene_stream)) == 0 )
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
    
    while ( bl_gff_read(&feature, gff_stream, BL_GFF_FIELD_ALL) == BL_READ_OK )
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
		int t;
		
		// Ensembl GFFs are not consistent with capitalization.
		// All lower case for Danio, capitalized for others
		//fprintf(stderr, "%d %zd %s\n", c, gene_names[c].count,
		//        BL_GFF_FEATURE_NAME(&feature));
		if ( (t = xt_alt_str_case_contains(&gene_names[c],
			BL_GFF_FEATURE_NAME(&feature))) >= 0 )
		{
		    // FIXME: Use accessor
		    printf("\n%s %s:\n", gff_basename, gene_names[c].strings[t]);
		    
		    // Path name format is important, parsed by other progs
		    // '-' is a separator, so don't allow it in species
		    // or gene names
		    strtr(gff_basename, "-", "_", 0);
		    
		    // Make gene name lower case in all files
		    // Don't enable until verifying that case in gene
		    // names has no meaning
		    // strlower(ALT_STR_STRINGS_AE(&gene_names[c], t));
		    strtr(ALT_STR_STRINGS_AE(&gene_names[c], t), "-", "_", 0);
		    snprintf(region_file, PATH_MAX,
			     "%s/%s-%s-%s-%" PRId64 "-%d-%" PRIu64 ".gff3",
			     output_dir, gff_basename,
			     ALT_STR_STRINGS_AE(&gene_names[c], t),
			     BL_GFF_SEQID(&feature), BL_GFF_START(&feature),
			     adjacent_genes, max_nt_distance);
		    status = extract_neighborhood(&feature, &gi, gff_stream,
						  header_stream, region_file,
						  adjacent_genes,
						  max_nt_distance);
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
 *      region_file.
 *
 *  Returns:
 *      EX_OK on success or another sysexits constant on error
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-13  Jason Bacon Begin
 ***************************************************************************/

int     extract_neighborhood(bl_gff_t *goi, bl_gff_index_t *gi,
	    FILE *gff_stream, FILE *header_stream, char *region_file,
	    int64_t adjacent_genes, int64_t max_nt_distance)

{
    bl_gff_t    neighbor;
    char        *neighbor_name;
    int64_t     g, len, distance;
    FILE        *region_stream;
    
    if ( (region_stream = fopen(region_file, "w")) == NULL )
    {
	fprintf(stderr, "ms-extract: Cannot open %s: %s\n",
		region_file, strerror(errno));
	return EX_CANTCREAT;
    }
    // Add the original GFF header to each neighborhood GFF
    bl_gff_copy_header(header_stream, region_stream);
    
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
    bl_gff_init(&neighbor);
    g = 0;
    // FIXME: Go adjacent_genes past the GOI in case the reverse seek
    // found fewer than adjacent_genes features within range
    while ( (g < adjacent_genes * 2 + 1) &&
	    (bl_gff_read(&neighbor, gff_stream, BL_GFF_FIELD_ALL) == BL_READ_OK) &&
	    ((strcmp(BL_GFF_TYPE(&neighbor), "###") == 0) ||
	     (strcmp(BL_GFF_SEQID(&neighbor), BL_GFF_SEQID(goi)) == 0)) )
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
	    
	    // FIXME: Restructure the loop so this is a loop condition
	    if ( distance > max_nt_distance )
		break;

	    len = BL_GFF_END(&neighbor) - BL_GFF_START(&neighbor);
	    
	    printf("%-3s %10" PRId64 " %10" PRId64 " %10" PRId64 " %10" PRId64 " %s\n",
		    BL_GFF_SEQID(&neighbor),
		    BL_GFF_START(&neighbor),
		    BL_GFF_END(&neighbor), len, distance, neighbor_name);
	
	    // Hack gene name into source field for DNA Feature Viewer.
	    // There should be a way to avoid this.
	    bl_gff_set_source_cpy(&neighbor, neighbor_name,
		BL_GFF_SOURCE_MAX_CHARS);
	    
	    // Output gene neighborhood in GFF format for above
	    bl_gff_write(&neighbor, region_stream,
		BL_GFF_FIELD_ALL);
	}
	bl_gff_free(&neighbor);
    }
    fclose(region_stream);
    
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
