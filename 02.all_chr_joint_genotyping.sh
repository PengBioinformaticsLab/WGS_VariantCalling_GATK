#!/bin/bash

WORKDIR=$1
sample_list=$2
chr_list=$3

cd $WORKDIR

for chr in $(cat $chr_list); do

  sh $WORKDIR/scripts/qsub_joint_genotyping.sh $WORKDIR $sample_list $chr ;

done
