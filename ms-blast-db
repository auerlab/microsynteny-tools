#!/bin/sh -e

# https://www.ncbi.nlm.nih.gov/books/NBK569850/?report=reader#!po=16.6667
# Database descriptions
# https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html
# landmark.tar.gz | Proteome of 27 model organisms
# nt.*tar.gz      | Partially non-redundant nucleotide sequences from 
#                   all traditional divisions of GenBank, EMBL, and DDBJ 
#                   excluding GSS,STS, PAT, EST, HTG, and WGS.
# landmark

# update_blastdb.pl --decompress nt

# makeblastdb -help|more

mkdir -p BLAST-DB
set -x
zcat Transcriptome/Homo_sapiens.GRCh38.cdna.all.fa.gz | \
    makeblastdb -title "Human Transcriptome" \
    -dbtype nucl -out BLAST-DB/Homo-sapiens-transcriptome

zcat Reference/Homo_sapiens.GRCh38.dna.toplevel.fa.xz | \
    makeblastdb -title "Human Genome" \
    -dbtype nucl -out BLAST-DB/Homo-sapiens-genome
