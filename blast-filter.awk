#############################################################################
#   Description:
#       Filter output of blastn -outfmt "6 qseqid pident evalue stitle"
#       to print "Query" "Identity" "E-value" "Location" "Gene".
#       For genome searches, there is no gene field, but it's OK
#       since $10 will just be blank.
#
#   History: 
#   Date        Name        Modification
#   2022-04-15  Jason Bacon Begin
#############################################################################

{
    split($7, a, ":");
    split($11, b, ":");
    chr="chr" a[3] ":" a[4] ":" a[5]
    printf("%-15s %10s %10s %10s %-25s %s\n", $1,$2,$3,$4,chr,b[2]);
}

