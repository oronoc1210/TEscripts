#!/bin/sh

export TEPATH="/global/projectb/scratch/cmodonog/TransposableElements"
export FILEPATH=$1
export FASTQ=$2
export OUTBASE=$3

echo "copying fastq file..."
time cp -v "$FILEPATH/$FASTQ" "$FILEPATH/$OUTBASE.fastq.gz"
echo "unzipping fastq file..."
time gunzip -v "$FILEPATH/$OUTBASE.fastq.gz"

time STAR --runMode alignReads --runThreadN 16 --winAnchorMultimapNmax 100 \
          --genomeDir "$TEPATH/Sbicolor_reference/STAR_index" \
          --readFilesIn "$FILEPATH/$OUTBASE.fastq" --outFilterMultimapNmax 100 \
          --outFileNamePrefix "$FILEPATH/$OUTBASE." --outSAMtype BAM SortedByCoordinate

rm "$FILEPATH/$OUTBASE.fastq"

time shifter --image=cmodonog/tetoolkit:v2.0.3 TEcount --project "$FILEPATH/$OUTBASE" \
             -b "$FILEPATH/$OUTBASE.Aligned.sortedByCoord.out.bam" --sortByPos --mode multi \
             --GTF "$TEPATH/Sbicolor_reference/Sbicolor_313_v3.1.gene.gtf" \
             --TE "$TEPATH/Sbicolor_reference/Sbicolor_313_v3.1.repeatmasked_assembly_v3.0.gtf" \
