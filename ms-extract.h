
#define MSYN_GENE_NAME_BUFF_LEN 256

#define DEFAULT_NT_DISTANCE     1000000
#define DEFAULT_ADJACENT_GENES  3

#define MS_OK                   0

void    usage(char *argv[]);
int     extract_neighborhood(bl_gff3_t *gene, bl_gff3_index_t *gi,
	    FILE *gff3_stream, FILE *header_stream, char *hood_file,
	    int64_t adjacent_genes, int64_t max_nt_distance);

