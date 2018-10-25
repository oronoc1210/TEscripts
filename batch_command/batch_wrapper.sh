#!/bin/sh

export TEPATH="/global/projectb/scratch/cmodonog/TransposableElements"
export FILEPATH=$1
export FASTQ=$2
export OUTBASE=$3

STAR --runMode alignReads --runThreadN 16 --winAnchorMultimapNmax 100 \
     --genomeDir "$TEPATH/Sbicolor_reference/STAR_index" \
     --readFilesIn "$TEPATH/$FILEPATH/$FASTQ" --outFilterMultimapNmax 100 \
     --outFileNamePrefix "$TEPATH/$FILEPATH/$OUTBASE" --outSAMtype BAM SortedByCoordinate

shifter image=cmodonog/tetoolkit:v2.0.3 TEcount --project "$TEPATH/$FILEPATH/$OUTBASE" \
    -b "$TEPATH/$FILEPATH/$OUTBASE.sortedByCoord.bam" --sortByPos --mode multi \
    --GTF "$TEPATH/Sbicolor_reference/Sbicolor_313_v3.1.gene.gtf" \
    --TE "$TEPATH/Sbicolor_reference/Sbicolor_313_v3.1.repeatmasked_assembly_v3.0.gtf" \
