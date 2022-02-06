/***************************************************************************
 *  Description:
 *      Show genes in the neighborhood of the specified gene for
 *      the specified species.
 *
 *  Arguments:
 *
 *  Returns:
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
#include "msyn.h"
#include "gff-index.h"

void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    bl_gff_t    gene = BL_GFF_INIT, neighbor_gene = BL_GFF_INIT;
    FILE        *gff_stream, *hood_stream, *header_stream;
    char        *gene_name,
		*gff_filename,
		*gff_basename,
		*end,
		*neighbor_name,
		hood_file[PATH_MAX],
		*output_dir = ".",
		*ext;
    uint64_t    max_nt_distance = DEFAULT_NT_DISTANCE,
		gene_count = DEFAULT_GENE_COUNT,
		g;
    int         arg;
    bl_gff_index_t    gi = BL_GFF_INDEX_INIT;

    if ( argc < 4 )
	    usage(argv);
    
    for (arg = 1; (arg < argc) && (argv[arg][0] == '-'); ++arg)
    {
	if ( strcmp(argv[arg], "--output-dir") == 0 )
	    output_dir = argv[++arg];
	else if ( strcmp(argv[arg], "--max-nt-distance") == 0 )
	{
	    max_nt_distance = strtoul(argv[++arg], &end, 10);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "msyn-hood: Invalid --max-nt-distance: %s\n",
			argv[arg]);
		usage(argv);
	    }
	}
	else if ( strcmp(argv[arg], "--gene-distance") == 0 )
	{
	    gene_count = strtoul(argv[++arg], &end, 10);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "msyn-hood: Invalid --gene-distance: %s\n",
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
    gene_name = argv[arg++];
    
    // FIXME: Maybe make xt_fopen_seekable() function to automatically
    // unzip compressed files if present
    if ( (gff_stream = fopen(gff_filename, "r")) == NULL )
    {
	fprintf(stderr, "msyn-hood: Could not open %s: %s\n",
		gff_filename, strerror(errno));
	exit(EX_NOINPUT);
    }
    header_stream = bl_gff_skip_header(gff_stream);
    while ( bl_gff_read(&gene, gff_stream, BL_GFF_FIELD_ALL) == BL_READ_OK )
    {
	if ( strcasecmp(BL_GFF_TYPE(&gene), "gene") == 0 )
	{
	    // Index positions of all genes in the file
	    if ( bl_gff_index_add_pos(&gi, BL_GFF_FILE_POS(&gene),
		 BL_GFF_SEQID(&gene), BL_GFF_START(&gene),
		 BL_GFF_END(&gene)) != BL_GFF_INDEX_OK )
	    {
		fprintf(stderr, "bl_gff_index_add_pos() failed.\n");
		return EX_SOFTWARE;
	    }
	    
	    if ( (BL_GFF_FEATURE_NAME(&gene) != NULL) &&
		 (strcasecmp(BL_GFF_FEATURE_NAME(&gene), gene_name) == 0) )
	    {
		printf("%s\t%" PRIu64 "\t%" PRIu64 "\t%s\n\n", BL_GFF_SEQID(&gene),
		    BL_GFF_START(&gene), BL_GFF_END(&gene), gene_name);
		
		if ( (gff_basename = strrchr(gff_filename, '/')) == NULL )
		    gff_basename = gff_filename;
		// Note: Destroys gff_filename == argv[1]
		if ( (ext = strchr(gff_basename, '.')) != NULL )
		    *ext = '\0';
		snprintf(hood_file, PATH_MAX, "%s/%s-%s.gff3",
		    output_dir, gff_basename, gene_name);
		if ( (hood_stream = fopen(hood_file, "w")) == NULL )
		{
		    fprintf(stderr, "msyn-hood: Cannot open %s: %s\n",
			    hood_file, strerror(errno));
		    return EX_CANTCREAT;
		}
		bl_gff_copy_header(header_stream, hood_stream);
		
		/*
		 *  Back up to at the position of this gene - max_nt_distance
		 *  and output all genes from there to position + max_nt_distance
		 */
		if ( bl_gff_index_seek_reverse(&gi, gff_stream, &gene,
			gene_count, max_nt_distance) != 0 )
		{
		    fprintf(stderr, "msyn-hood: Seek %zu failed.\n", BL_GFF_FILE_POS(&gene));
		    return EX_SOFTWARE;
		}
		
		// From leftmost neighbor read gene_count before and after ref
		for (g = 0; (g < gene_count * 2 + 1) &&
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
				BL_GFF_END(&gene) + max_nt_distance )
			    break;
		    }
		    bl_gff_free(&neighbor_gene);
		}
		fclose(hood_stream);
		
		// FIXME: Temporary hack
		break;
	    }
	}
	bl_gff_free(&gene);
    }
    fclose(gff_stream);
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s\n"
		    "   [--output-dir dir]\n"
		    "   [--gene-distance genes]\n"
		    "   [--max-nt-distance nucleotides]\n"
		    "   file.gff gene-name\n", argv[0]);
    exit(EX_USAGE);
}

