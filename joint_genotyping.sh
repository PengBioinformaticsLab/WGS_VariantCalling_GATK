#!/bin/bash
#SBATCH -A r00346
#SBATCH --mail-user=jiaji@iu.edu
#SBATCH --array=1-24  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH -J GenomicsDBImport
#SBATCH -o logs/GenomicsDBImport_chr%a.out

WORKDIR=$1
chr_list=$2
sample_gvcf_map=$3
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
mapfile -t chromosomes < $chr_list
chr=${chromosomes[$SLURM_ARRAY_TASK_ID - 1]}

[[ -f "${sample_gvcf_map}" ]] || die "Sample map file not found: ${sample_gvcf_map}"
[[ -n "${chr}" ]] || die "Chromosome not found in list"

set -eo pipefail

# Step 1: GenomicsDBImport
if [ ! -d "$WORKDIR/GVCF_DB/${chr}_DB" ]; then
	echo_step "Aggregating gVCFs for $chr"
	gatk --java-options "-Xmx32g -XX:+UseParallelGC -XX:ParallelGCThreads=12 -XX:GCTimeLimit=50" \
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
	gatk --java-options "-Xmx32g -XX:+UseParallelGC -XX:ParallelGCThreads=8" \
		GenotypeGVCFs \
		-R $REF \
		-V gendb://$WORKDIR/GVCF_DB/${chr}_DB \
		-L $chr \
		-O $WORKDIR/variant_calls/vcf/${chr}.vcf.gz \
		--tmp-dir $WORKDIR/temp || die "GenotypeGVCFs failed for ${chr}"
fi

echo_step "Completed GenomicsDBImport and GenotypeGVCFs for ${chr}"