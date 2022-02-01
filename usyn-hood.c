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
#include "findex.h"

void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    bl_gff_t    feature = BL_GFF_INIT;
    FILE        *gff_stream;
    char        *gene_name,
		*gff_filename,
		gene_name_string[USYN_GENE_NAME_BUFF_LEN],
		*end;
    size_t      distance;
    findex_t    findex = FINDEX_INIT;
    
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
    findex_add_pos(&findex, gff_stream);
    while ( bl_gff_read(&feature, gff_stream, BL_GFF_FIELD_ALL) == BL_READ_OK )
    {
	if ( (strcmp(BL_GFF_FEATURE(&feature), "gene") == 0) &&
	     (strstr(BL_GFF_ATTRIBUTES(&feature), gene_name_string) != NULL) )
	{
	    bl_gff_write(&feature, stdout, BL_GFF_FIELD_ALL);
	    
	    /*
	     *  Back up to at the position of this gene - distance
	     *  and output all genes from there to position + distance
	     */
	    
	}
    }
    fclose(gff_stream);
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: file.gff gene-name distance%s\n", argv[0]);
    exit(EX_USAGE);
}
