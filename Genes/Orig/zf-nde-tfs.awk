#############################################################################
#   Description:
#  
#   Arguments:
#
#   Returns:
#
#   History: 
#   Date        Name        Modification
#   2022-03-23  Jason Bacon Begin
#############################################################################

BEGIN {
    IFS="[\t]";
}
$1 != "ext_gene" {
    split($1, a, " ");
    if ( ($1 != $2) && ($2 != "N/A") )
    {
	printf("%s|%s\n", a[1], $2);
    }
    else
    {
	print a[1];
    }
}
