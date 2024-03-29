
#ifndef _BIOLIBC_GFF3_H_
#include <biolibc/gff3.h>
#endif

#ifndef __bool_true_false_are_defined
#include <stdbool.h>
#endif

typedef struct
{
    size_t      array_size;
    size_t      count;
    size_t      goi_index;
    size_t      goi_len;
    bl_gff3_t    *features;  // Array of GFF features
    char        *species;   // FIXME: Maybe use fixed size strings
    char        *goi;
    char        *chrom;
}   bl_gff3_region_t;

#include "gff-region-rvs.h"
#include "gff-region-accessors.h"
#include "gff-region-mutators.h"

/* gff-region.c */
int bl_gff3_region_load(bl_gff3_region_t *region, const char *filename);
void bl_gff3_region_init(bl_gff3_region_t *region);
void bl_gff3_region_free(bl_gff3_region_t *region);
int bl_gff3_region_commonality(bl_gff3_region_t *h1, bl_gff3_region_t *h2);
bl_gff3_region_t *bl_gff3_region_intersect(bl_gff3_region_t *r1, bl_gff3_region_t *r2);
bool bl_gff3_region_duplicate_gene(bl_gff3_region_t *r, const char *feature_name);

