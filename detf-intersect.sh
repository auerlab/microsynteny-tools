#!/bin/sh -e

for gene in $(cat Genes/detf-names.txt); do
    ./intersect.sh "$gene"
done
