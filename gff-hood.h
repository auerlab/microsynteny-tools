#define BL_GFF_HOOD_MAX_FEATURES   64

typedef struct
{
    size_t      gene_count;
    size_t      goi_index;
    bl_gff_t    features[BL_GFF_HOOD_MAX_FEATURES];
    char        *species;
    char        *goi;
}   bl_gff_hood_t;

/* gff-hood.c */
int bl_gff_hood_load(bl_gff_hood_t *hood, const char *filename);
void bl_gff_hood_init(bl_gff_hood_t *hood);
void bl_gff_hood_free(bl_gff_hood_t *hood);
int bl_gff_hood_commonality(bl_gff_hood_t *h1, bl_gff_hood_t *h2);
