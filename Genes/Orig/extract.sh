#!/bin/sh -e

awk '{ print $2 }' detf-refined.tsv > ../detf-refined-names.txt

awk -F , 'length($2) > 0 && $1 ~ "ENSDART" { print $2 }' all-de.csv \
    | awk -f split-gene.awk > ../all-de-names.txt

