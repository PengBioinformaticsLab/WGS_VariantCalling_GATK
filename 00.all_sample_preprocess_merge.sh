#!/bin/bash

# Define the directory containing BAM files and working directory
 # Update with the correct paths
WORKDIR=$1
sample_list=$2

# Navigate to the working directory
cd $WORKDIR

# Loop through each BAM file in the directory
for SAMPLE in $(cat $sample_list); do

  sh $WORKDIR/scripts/qsub_preprocess_merge.sh $WORKDIR $SAMPLE ;

done

## copy the unique BAMs to the mergedBAMs directory
#for SAMPLE in $(cat BAM/samples.uniq); do
# cp BAM/*/${SAMPLE}.hg38_bwa.sorted.bam mergedBAMs/${SAMPLE}.sorted.bam
# cp BAM/*/${SAMPLE}.hg38_bwa.sorted.bam.bai mergedBAMs/${SAMPLE}.sorted.bam.bai
#done