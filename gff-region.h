typedef struct
{
    size_t      count;
    size_t      goi_index;
    bl_gff_t    *features;  // Array of GFF features
    char        *species;
    char        *goi;
}   bl_gff_region_t;

/* gff-region.c */
int bl_gff_region_load(bl_gff_region_t *region, const char *filename);
void bl_gff_region_init(bl_gff_region_t *region);
void bl_gff_region_free(bl_gff_region_t *region);
int bl_gff_region_commonality(bl_gff_region_t *h1, bl_gff_region_t *h2);
bl_gff_region_t *bl_gff_region_intersect(bl_gff_region_t *r1, bl_gff_region_t *r2);

