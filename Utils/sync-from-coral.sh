#!/bin/sh -e

for dir in GFF Transcriptome Reference BLAST-DB; do
    rsync -av --partial --progress ~/Coral/Prog/Src/microsynteny-tools/$dir .
done
