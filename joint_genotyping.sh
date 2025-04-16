#!/bin/bash

WORKDIR=$1
sample_gvcf_map=$2
chr_list=$3

###### match the WGS sample names with the RNAseq sample names, set col1 as sample names in the $sample_gvcf_map

module load gatk
module load java/17.0.7

REF=/N/slate/jiaji/CMG_references/hg38/GATK4_bundle_resource/Homo_sapiens_assembly38.fasta

# --------------------------
# Functions
# --------------------------

set -e

function die() {
    echo "ERROR: $*" >&2
    exit 1
}

function echo_step() {
	echo "########################################"
    echo -e "# $(date '+%m-%d-%Y %H:%M:%S') - $1 #"
    echo "########################################"
}

cd $WORKDIR
mkdir -p $WORKDIR/temp $WORKDIR/GVCF_DB $WORKDIR/variant_calls/vcf
rm -rf $WORKDIR/temp/*

# --------------------------
# Chromosome Processing
# --------------------------

set -eo pipefail

# Step 1: GenomicsDBImport
if [ ! -d "$WORKDIR/GVCF_DB/${chr}_DB" ]; then
	echo_step "Aggregating gVCFs for $chr"
	gatk --java-options "-Xmx32g -XX:+UseParallelGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50" \
		GenomicsDBImport \
		--genomicsdb-workspace-path $WORKDIR/GVCF_DB/${chr}_DB \
		--sample-name-map $sample_gvcf_map \
		-L $chr \
		--reader-threads 4 \
		--tmp-dir $WORKDIR/temp || die "GenomicsDBImport failed for ${chr}"
else
  echo_step "Existing database found for ${chr}, skipping import"
fi


# Step 2: GenotypeGVCFs
if [[ ! -s $WORKDIR/variant_calls/vcf/${chr}.vcf.gz ]]; then
	echo_step "Joint Genotyping for $chr"
	gatk --java-options "-Xmx32g -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
		GenotypeGVCFs \
		-R $REF \
		-V gendb://$WORKDIR/GVCF_DB/${chr}_DB \
		-L $chr \
		-O $WORKDIR/variant_calls/vcf/${chr}.vcf.gz \
		--tmp-dir $WORKDIR/temp || die "GenotypeGVCFs failed for ${chr}"
fi

echo_step "Completed GenomicsDBImport and GenotypeGVCFs for ${chr}"