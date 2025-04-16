#!/bin/bash

WORKDIR=$1
sample_list=$2

# Navigate to the working directory
cd $WORKDIR

# Loop through each BAM file in the directory
for SAMPLE in $(cat $sample_list); do

  sh $WORKDIR/scripts/qsub_markdup_variantCall.sh $WORKDIR $SAMPLE ;

done
