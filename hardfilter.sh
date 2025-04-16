#!/bin/bash

#SBATCH -A r00302
#SBATCH --mail-user=jiaji@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --mem=4gb
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH -J WGS_gatk_filtering
#SBATCH -o log/WGS_gatk_filtering_%j.out

WORKDIR=$1

cd $WORKDIR

module load gatk
module load bcftools
module load java/17.0.7
module load vcftools

picard=/N/soft/rhel8/picard/3.0.0/build/libs/picard.jar
export PATH="/N/slate/jiaji/tools/annovar:$PATH"

WD=${WORKDIR}/WGS_GATK/filtration


export OMP_NUM_THREADS=4

logfile="$WORKDIR/log/WGS_GATK_$(date +"%Y%m%d%H%M").log"

set -ex
exec 2>>$logfile



start_time_mod=$(date +%s)

if [ ! -s $WD/WGS.GATK.joint.vqsr.snps.vcf.gz ];then
	gatk SelectVariants \
		-V $WD/WGS.GATK.joint.vqsr.vcf.gz \
		-select-type SNP \
		-O $WD/WGS.GATK.joint.vqsr.snps.vcf.gz
fi

if [ ! -s $WD/WGS.GATK.joint.vqsr.indels.vcf.gz ];then
	gatk SelectVariants \
		-V $WD/WGS.GATK.joint.vqsr.vcf.gz \
		-select-type INDEL \
		-O $WD/WGS.GATK.joint.vqsr.indels.vcf.gz
fi

if [ ! -s $WD/WGS.GATK.joint.snps.filtered.vcf.gz ];then
	gatk VariantFiltration \
		-V $WD/WGS.GATK.joint.vqsr.snps.vcf.gz \
		-filter "DP < 4" --filter-name "DP4" \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "SOR > 3.0" --filter-name "SOR3" \
		-filter "FS > 60.0" --filter-name "FS60" \
		-filter "MQ < 40.0" --filter-name "MQ40" \
		-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
		-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
		-O $WD/WGS.GATK.joint.snps.filtered.vcf.gz
fi

if [ ! -s $WD/WGS.GATK.joint.indels.filtered.vcf.gz ]; then
	gatk VariantFiltration \
		-V $WD/WGS.GATK.joint.vqsr.indels.vcf.gz \
		-filter "DP < 4" --filter-name "DP4" \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "FS > 200.0" --filter-name "FS200" \
		-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-8" \
		-O $WD/WGS.GATK.joint.indels.filtered.vcf.gz
fi


if [ ! -s $WD/WGS.GATK.joint.vqsr.filtered.snps.pass.vcf.gz.tbi ]; then
	vcftools --gzvcf $WD/WGS.GATK.joint.snps.filtered.vcf.gz \
	--remove-filtered-all \
	--recode \
	--stdout | bgzip -@ 4 -c > $WD/WGS.GATK.joint.vqsr.filtered.snps.pass.vcf.gz
	tabix -p vcf $WD/WGS.GATK.joint.vqsr.filtered.snps.pass.vcf.gz
fi

if [ ! -s $WD/WGS.GATK.joint.vqsr.filtered.indels.pass.vcf.gz.tbi ]; then
	vcftools --gzvcf $WD/WGS.GATK.joint.indels.filtered.vcf.gz \
	--remove-filtered-all \
	--recode \
	--stdout | bgzip -@ 4 -c > $WD/WGS.GATK.joint.vqsr.filtered.indels.pass.vcf.gz
	tabix -p vcf $WD/WGS.GATK.joint.vqsr.filtered.indels.pass.vcf.gz
fi

end_time_mod=$(date +%s)
if [ "$OSTYPE" = "darwin"* ]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [ "$OSTYPE" = "darwin"* ]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Hard Filtration Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec">>$logfile

