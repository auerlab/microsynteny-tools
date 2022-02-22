#!/bin/sh -e

for gene in $(cat DETF-names.txt); do
    ./divergence.sh $gene
done
