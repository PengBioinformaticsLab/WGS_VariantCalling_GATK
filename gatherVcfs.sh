#!/bin/bash

#SBATCH -A r00302
#SBATCH --mail-user=jiaji@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=90:00:00
#SBATCH --mem=64gb
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH -J WGS_MergeVcfs
#SBATCH -o log/WGS_gatk_mergeVcfs_%j.out

export OMP_NUM_THREADS=8

set -euo pipefail

module load gatk
module load java/17.0.7

BUNDLE_RESOURCE=/N/slate/jiaji/CMG_references/hg38/GATK4_bundle_resource
REF=$BUNDLE_RESOURCE/Homo_sapiens_assembly38.fasta
DBSNP=$BUNDLE_RESOURCE/Homo_sapiens_assembly38.dbsnp138.vcf

WORKDIR=$1
cd $WORKDIR

gatk --java-options "-Xmx32g" MergeVcfs \
    -I vcf_list.txt \
    -O $WORKDIR/variant_calls/WGS_GATK.merged.vcf.gz

# Final Validation and Indexing 
gatk --java-options "-Xmx32g" IndexFeatureFile \
    -I $WORKDIR/variant_calls/WGS_GATK.merged.vcf.gz \
	--num-threads 8

gatk --java-options "-Xmx32g" ValidateVariants \
    -V $WORKDIR/variant_calls/WGS_GATK.merged.vcf.gz \
    -R $REF \
    --dbsnp $DBSNP
