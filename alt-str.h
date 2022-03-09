
/*
 *  List of equivalent genes.  genes[0] is the primary,
 *  genes[1] secondary, etc.  I.e. if genes[0] is not present,
 *  genes[1] can be used.
 */

typedef struct
{
    size_t  count;
    char    **strings;
}   alt_str_t;

/* alt-str.c */                                                                 
int xt_alt_str_read_line_malloc(alt_str_t *alt_str, FILE *stream);
ssize_t xt_alt_str_inhale_list(alt_str_t **list, FILE *stream);
int xt_alt_str_case_contains(alt_str_t *alt_str, char *str);
