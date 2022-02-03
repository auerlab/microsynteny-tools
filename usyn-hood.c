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
#include <biolibc/gff.h>
#include "usyn.h"
#include "gff-index.h"

void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    bl_gff_t    gene = BL_GFF_INIT, neighbor_gene = BL_GFF_INIT;
    FILE        *gff_stream;
    char        *gene_name,
		*gff_filename,
		gene_name_string[USYN_GENE_NAME_BUFF_LEN],
		*end,
		*neighbor_name;
    size_t      distance;
    bl_gff_index_t    gi = BL_GFF_INDEX_INIT;
    
    switch(argc)
    {
	case 4:
	    break;
	
	default:
	    usage(argv);
    }
    gff_filename = argv[1];
    gene_name = argv[2];
    
    snprintf(gene_name_string, USYN_GENE_NAME_BUFF_LEN, "Name=%s;", gene_name);
    distance = strtoul(argv[3], &end, 10);
    if ( *end != '\0' )
    {
	fprintf(stderr, "usyn-hood: Invalid distance: %s\n", argv[3]);
	usage(argv);
    }

    // FIXME: Maybe make xt_fopen_seekable() function to automatically
    // unzip compressed files if present
    if ( (gff_stream = fopen(gff_filename, "r")) == NULL )
    {
	fprintf(stderr, "usyn-hood: Could not open %s: %s\n",
		gff_filename, strerror(errno));
	exit(EX_NOINPUT);
    }
    bl_gff_skip_header(gff_stream);
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
		
		/*
		 *  Back up to at the position of this gene - distance
		 *  and output all genes from there to position + distance
		 */
		if ( bl_gff_index_seek_first_ge(&gi, gff_stream,
			BL_GFF_SEQUENCE(&gene),
			BL_GFF_END(&gene) - distance) != 0 )
		{
		    fprintf(stderr, "usyn-hood: Seek %zu failed.\n", BL_GFF_FILE_POS(&gene));
		    return EX_SOFTWARE;
		}
		while ( bl_gff_read(&neighbor_gene, gff_stream,
				    BL_GFF_FIELD_ALL) == BL_READ_OK )
		{
		    //printf("%" PRIu64 "\n", BL_GFF_START(&neighbor_gene));
		    if ( strcmp(BL_GFF_FEATURE(&neighbor_gene), "gene") == 0 )
		    {
			neighbor_name =
			    strstr(BL_GFF_ATTRIBUTES(&neighbor_gene), "Name=");
			if ( neighbor_name != NULL )
			{
			    end = strchr(neighbor_name, ';');
			    *end = '\0';
			}
			else
			    neighbor_name = "unnamed";
			printf("%s\t%" PRIu64 "\t%" PRIu64 "\t%s\n",
				BL_GFF_SEQUENCE(&neighbor_gene),
				BL_GFF_START(&neighbor_gene),
				BL_GFF_END(&neighbor_gene), neighbor_name);
			//printf("%" PRIu64 " %" PRIu64 "\n", BL_GFF_START(&neighbor_gene),
			//    BL_GFF_START(&gene) + distance);
			if ( BL_GFF_START(&neighbor_gene) > BL_GFF_START(&gene) + distance )
			    break;
		    }
		    bl_gff_free(&neighbor_gene);
		}
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
