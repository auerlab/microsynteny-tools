#############################################################################
#   Description:
#  
#   Arguments:
#
#   Returns:
#
#   History: 
#   Date        Name        Modification
#   2022-03-30  Jason Bacon Begin
#############################################################################

BEGIN {
    
}
{
    genes += $1;
    conserved += $2;
    changed += $3;
}
END {
    printf("Genes: %d  Conserved: %d (%0.1f%%)  Changed: %d (%0.1f%%)\n",
	    genes, conserved, conserved / genes * 100,
	    changed, changed / genes * 100);
}

