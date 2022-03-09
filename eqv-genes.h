
/*
 *  List of equivalent genes.  genes[0] is the primary,
 *  genes[1] secondary, etc.  I.e. if genes[0] is not present,
 *  genes[1] can be used.
 */

typedef struct
{
    size_t  count;
    char    **genes;
}   eqv_genes_t;
