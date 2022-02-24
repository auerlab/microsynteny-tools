#!/bin/sh -e

for gene in $(cat DETF-names.txt); do
    ./common-to.sh $gene
done
