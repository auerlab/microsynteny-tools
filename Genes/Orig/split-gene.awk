#############################################################################
#   Description:
#  
#   Arguments:
#
#   Returns:
#
#   History: 
#   Date        Name        Modification
#   2022-02-28  Jason Bacon Begin
#############################################################################

{
    # Assuming / means "or" in gene names, e.g.
    # il17a/f3 means il17a or il17f3
    if ( $1 ~ "/" ) {
	split($1, a, "/");
	print a[1]; print substr(a[1], 1, length(a[1]) - 1) a[2];
    } else {
	print $1;
    }
}
