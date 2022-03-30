#!/bin/sh -e

# https://www.ncbi.nlm.nih.gov/books/NBK569850/?report=reader#!po=16.6667
# Database descriptions
# https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html
# landmark.tar.gz | Proteome of 27 model organisms
# nt.*tar.gz      | Partially non-redundant nucleotide sequences from 
#                   all traditional divisions of GenBank, EMBL, and DDBJ 
#                   excluding GSS,STS, PAT, EST, HTG, and WGS.
update_blastdb.pl --decompress nt landmark
