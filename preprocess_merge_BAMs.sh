#!/bin/bash

BAM1=$1
BAM2=$2
SAMPLE=$3
WORKDIR=$4

module load gatk
module load java/17.0.7
module load samtools

picard=/N/soft/rhel8/picard/3.0.0/build/libs/picard.jar


if [ ! -s $WORKDIR/mergedBAMs/${SAMPLE}.sorted.bam.bai ];then
    java -Xmx24G -jar $picard MergeSamFiles \
        -I $BAM1 -I $BAM2 \
		-O $WORKDIR/mergedBAMs/${SAMPLE}.sorted.bam USE_THREADING=true --CREATE_INDEX true
fi

