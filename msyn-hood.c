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
		gene_name_string[MSYN_GENE_NAME_BUFF_LEN],
		*end,
		*neighbor_name,
		hood_file[PATH_MAX],
		*output_dir = ".";
    size_t      distance;
    int         arg;
    bl_gff_index_t    gi = BL_GFF_INDEX_INIT;

    if ( argc < 4 )
	    usage(argv);
    
    for (arg = 1; (arg < argc) && (argv[arg][0] == '-'); ++arg)
    {
	if ( strcmp(argv[arg], "--output-dir") == 0 )
	    output_dir = argv[++arg];
    }
    
    if ( arg > argc - 3 )
	usage(argv);
    gff_filename = argv[arg++];
    gene_name = argv[arg++];
    
    snprintf(gene_name_string, MSYN_GENE_NAME_BUFF_LEN, "Name=%s;", gene_name);
    distance = strtoul(argv[arg++], &end, 10);
    if ( *end != '\0' )
    {
	fprintf(stderr, "msyn-hood: Invalid distance: %s\n", argv[3]);
	usage(argv);
    }

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
	if ( strcmp(BL_GFF_FEATURE(&gene), "gene") == 0 )
	{
	    // Index positions of all genes in the file
	    bl_gff_index_add_pos(&gi, BL_GFF_FILE_POS(&gene),
		BL_GFF_SEQUENCE(&gene), BL_GFF_START(&gene), BL_GFF_END(&gene));
	    
	    if ( strstr(BL_GFF_ATTRIBUTES(&gene), gene_name_string) != NULL )
	    {
		printf("%s\t%" PRIu64 "\t%" PRIu64 "\t%s\n\n", BL_GFF_SEQUENCE(&gene),
		    BL_GFF_START(&gene), BL_GFF_END(&gene), gene_name);
		
		if ( (gff_basename = strrchr(gff_filename, '/')) == NULL )
		    gff_basename = gff_filename;
		snprintf(hood_file, PATH_MAX, "%s/%s-%s-hood.gff3",
		    output_dir, gff_basename, gene_name);
		if ( (hood_stream = fopen(hood_file, "w")) == NULL )
		{
		    fprintf(stderr, "msyn-hood: Cannot open %s: %s\n",
			    hood_file, strerror(errno));
		    return EX_CANTCREAT;
		}
		bl_gff_copy_header(header_stream, hood_stream);
		
		/*
		 *  Back up to at the position of this gene - distance
		 *  and output all genes from there to position + distance
		 */
		if ( bl_gff_index_seek_first_ge(&gi, gff_stream,
			BL_GFF_SEQUENCE(&gene),
			BL_GFF_END(&gene) - distance) != 0 )
		{
		    fprintf(stderr, "msyn-hood: Seek %zu failed.\n", BL_GFF_FILE_POS(&gene));
		    return EX_SOFTWARE;
		}
		while ( bl_gff_read(&neighbor_gene, gff_stream,
				    BL_GFF_FIELD_ALL) == BL_READ_OK )
		{
		    if ( strcmp(BL_GFF_FEATURE(&neighbor_gene), "gene") == 0 )
		    {
			neighbor_name =
			    strstr(BL_GFF_ATTRIBUTES(&neighbor_gene), "Name=");
			if ( neighbor_name != NULL )
			{
			    neighbor_name = strchr(neighbor_name, '=') + 1;
			    end = strchr(neighbor_name, ';');
			    *end = '\0';
			}
			else
			    neighbor_name = "unnamed";
			printf("%s\t%" PRIu64 "\t%" PRIu64 "\t%s\n",
				BL_GFF_SEQUENCE(&neighbor_gene),
				BL_GFF_START(&neighbor_gene),
				BL_GFF_END(&neighbor_gene), neighbor_name);
		    
			// Hack gene name into source field for
			// DNA Feature Viewer.  Should be a way to avoid this.
			bl_gff_set_source_cpy(&neighbor_gene, neighbor_name,
			    BL_GFF_SOURCE_MAX_CHARS);
			// Output gene neighborhood in GFF format for above
			bl_gff_write(&neighbor_gene, hood_stream,
			    BL_GFF_FIELD_ALL);
			
			if ( BL_GFF_START(&neighbor_gene) > BL_GFF_START(&gene) + distance )
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
    fprintf(stderr, "Usage: file.gff gene-name distance%s\n", argv[0]);
    exit(EX_USAGE);
}

