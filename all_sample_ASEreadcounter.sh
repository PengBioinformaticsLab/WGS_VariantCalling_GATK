#!/bin/bash

WORKDIR=$1
sample_bam_map=$2
VCF=$3

cd $WORKDIR


while IFS=$'\t' read -r SAMPLE BAM; do 
    sh scripts/qsub_ASEReadCounter.sh $WORKDIR $SAMPLE $BAM $VCF ; 
done < $sample_bam_map

