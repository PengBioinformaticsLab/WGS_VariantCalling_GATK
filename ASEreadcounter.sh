#!/bin/bash

WORKDIR=$1
SAMPLE=$2
BAM=$3
VCF=$4

cd $WORKDIR

module load gatk
module load java/17.0.7
module load bcftools

REF=/N/slate/jiaji/CMG_references/hg38/GATK4_bundle_resource/Homo_sapiens_assembly38.fasta
WD=$WORKDIR/ASE/${SAMPLE}
mkdir -p $WD

logfile="$WORKDIR/log/${SAMPLE}.ASE.$(date +"%Y%m%d%H%M").log"

set -ex
exec 2>>$logfile

function echo_step() {
    echo "########################################"
    echo -e "# $(date '+%m-%d-%Y %H:%M:%S') - $1 "
    echo "########################################"
}

######## Optional, match RG in the bam files with column names in the joint vcf
if [ ! -s $WD/${SAMPLE}.star.sorted.bam ]; then
	gatk AddOrReplaceReadGroups \
		-I $BAM \
		-O $WD/${SAMPLE}.star.sorted.bam \
		-RGID ${SAMPLE} -RGLB mRNAseq -RGPL ILLUMINA -RGSM ${SAMPLE} -RGPU AHW3MGDMXX --CREATE_INDEX
fi

# keep only the heterozygous sites
if [ ! -s $WD/${SAMPLE}.gatk.filtered.pass.snps.het.vcf.gz ]; then
	echo_step "BCFtools Filtering... "
	bcftools view -s $SAMPLE --threads 4 $VCF -Ou | bcftools view -i 'GT="0/1"' --threads 4 -Oz -o $WD/${SAMPLE}.gatk.filtered.pass.snps.het.vcf.gz
	tabix -p vcf $WD/${SAMPLE}.gatk.filtered.pass.snps.het.vcf.gz
fi

# Alele-specific expression (ASE) read counter
# Counts the number of reads supporting each allele at a given site in a BAM file
# Only counts on the Biallelic sites
if [ ! -s $WD/${SAMPLE}.ase.table ]; then
	echo_step "Generating ASE table"
	gatk --java-options "-Xmx32g" ASEReadCounter \
		-R $REF \
		-I $WD/${SAMPLE}.star.sorted.bam \
		-V $WD/${SAMPLE}.gatk.filtered.pass.snps.het.vcf.gz \
		-O $WD/${SAMPLE}.ase.table
fi


