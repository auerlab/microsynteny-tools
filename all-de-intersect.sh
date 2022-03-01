#!/bin/sh -e

for gene in $(cat Genes/all-de-names.txt); do
    ./intersect.sh $gene
done
