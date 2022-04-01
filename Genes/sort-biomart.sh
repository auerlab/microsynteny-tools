#!/bin/sh -e

raw=biomart-orthologs.tsv.xz
sorted=biomart-orthologs-sorted.tsv
xzcat $raw | grep -m 1 '^Gene name' > $sorted
xzcat $raw | grep -v '^Gene name' | sort --ignore-case | uniq > $sorted
