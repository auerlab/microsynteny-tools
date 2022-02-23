typedef struct
{
    size_t      count;
    size_t      goi_index;
    bl_gff_t    *features;  // Array of GFF features
    char        *species;
    char        *goi;
}   bl_gff_region_t;

/* gff-hood.c */
int bl_gff_region_load(bl_gff_region_t *hood, const char *filename);
void bl_gff_region_init(bl_gff_region_t *hood);
void bl_gff_region_destroy(bl_gff_region_t *hood);
int bl_gff_region_commonality(bl_gff_region_t *h1, bl_gff_region_t *h2);
bl_gff_region_t *bl_gff_region_intersect(bl_gff_region_t *r1, bl_gff_region_t *r2);

