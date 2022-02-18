#!/bin/sh -e

for gene in $(awk '{ print $2 }' DETF_refined.tsv); do
    printf "===\n$gene\n"
    ./divergence.sh $gene
done
