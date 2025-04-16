#!/bin/bash

# --------------------------
# Input Parameters
# --------------------------
WORKDIR=$1
SAMPLE=$2

# --------------------------
# Module Load & Environment
# --------------------------
module load gatk
module load java/17.0.7
#module load conda
#conda activate r_env

# --------------------------
# Path Configuration
# --------------------------
picard=/N/soft/rhel8/picard/3.0.0/build/libs/picard.jar
BUNDLE_RESOURCE=/N/slate/jiaji/CMG_references/hg38/GATK4_bundle_resource
REF=$BUNDLE_RESOURCE/Homo_sapiens_assembly38.fasta
DBSNP=$BUNDLE_RESOURCE/Homo_sapiens_assembly38.dbsnp138.vcf
KNOWNINDEL=$BUNDLE_RESOURCE/Homo_sapiens_assembly38.known_indels.vcf.gz
MILLS=$BUNDLE_RESOURCE/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# --------------------------
# Functions
# --------------------------
function die() {
    echo "ERROR: $*" >&2
    exit 1
}

function echo_step() {
	echo "########################################"
    echo -e "# $(date '+%m-%d-%Y %H:%M:%S') - $1 #"
    echo "########################################"
}

# --------------------------
# Directory Setup
# --------------------------
mkdir -p "${WORKDIR}/output/bqsr" \
         "${WORKDIR}/table" \
         "${WORKDIR}/variant_calls/gvcf" 


# --------------------------
# Pipeline Steps
# --------------------------
set -eo pipefail

# Step 1: Mark Duplicates 3hrs
if [[ ! -s $WORKDIR/table/${SAMPLE}.marked_dup_metrics.txt ]]; then
    echo_step "Mark Duplicates..."
	gatk --java-options "-Xmx16g" \
		MarkDuplicatesSpark \
		-I $WORKDIR/mergedBAMs/${SAMPLE}.sorted.bam \
		-O $WORKDIR/output/${SAMPLE}.sort.dup.bam \
		-M $WORKDIR/table/${SAMPLE}.marked_dup_metrics.txt  \
		--ASSUME_SORT_ORDER coordinate \
		--create-output-bam-index true || die "MarkDuplicates failed"
fi

# Step 2: Base Recalibration 3-5hrs
if [[ ! -s -O $WORKDIR/table/${SAMPLE}.recalBefore.table]]; then
	echo_step "Calculating BQSR Table..."
	gatk --java-options "-Xmx16g" \
		BaseRecalibrator \
		-I $WORKDIR/output/${SAMPLE}.sort.dup.bam \
		-R $REF \
		--known-sites $DBSNP \
		--known-sites $KNOWNINDEL \
		--known-sites $MILLS \
		-O $WORKDIR/table/${SAMPLE}.recalBefore.table || die "BaseRecalibrator failed"
fi

# Step 3: Apply BQSR 3-4hrs
if [[ ! -s $WORKDIR/output/bqsr/${SAMPLE}.sort.dup.bqsr.bam ]]; then
	echo_step "Applying Recalibration..."
	gatk --java-options "-Xmx16g" \
		ApplyBQSR \
		-R $REF \
		-I $WORKDIR/output/${SAMPLE}.sort.dup.bam \
		--bqsr-recal-file $WORKDIR/table/${SAMPLE}.recalBefore.table \
		-O $WORKDIR/output/bqsr/${SAMPLE}.sort.dup.bqsr.bam || die "ApplyBQSR failed"
fi

# Step 4: Variant Calling
# 20hrs
if [[ ! -s $WORKDIR/variant_calls/gvcf/${SAMPLE}.g.vcf.gz ]]; then
	echo_step "Calling Variants..."
	gatk --java-options "-Xmx16g" \
		HaplotypeCaller \
		-I $WORKDIR/output/bqsr/${SAMPLE}.sort.dup.bqsr.bam \
		-R $REF \
		-ERC GVCF \
		-O $WORKDIR/variant_calls/gvcf/${SAMPLE}.g.vcf.gz  || die "HaplotypeCaller failed"
fi


echo_step "Completed"
