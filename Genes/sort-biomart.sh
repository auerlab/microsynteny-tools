#!/bin/sh -e

raw=biomart-orthologs.tsv.xz
sorted=biomart-orthologs-sorted.tsv

# Copy header
xzcat $raw | grep -m 1 '^Gene name' > $sorted

# Sort all except header
xzcat $raw | grep -v '^Gene name' | sort --ignore-case | uniq | xz > $sorted

xz $sorted
