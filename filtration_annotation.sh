#!/bin/bash

#SBATCH -A r00302
#SBATCH --mail-user=jiaji@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=96:00:00
#SBATCH --mem=64gb
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH -J WGS_gatk_filtering
#SBATCH -o log/WGS_gatk_filtering_%j.out

WORKDIR=$1
JOINTVCF=$2

cd $WORKDIR

module load gatk
module load bcftools
module load java/17.0.7
module load vcftools

picard=/N/soft/rhel8/picard/3.0.0/build/libs/picard.jar
export PATH="/N/slate/jiaji/tools/annovar:$PATH"

WD=${WORKDIR}/WGS_GATK/filtration
WD_anno=$(dirname "$WD")/annovar

mkdir -p $WD $WD_anno

## generate log  file 
logfile="$WORKDIR/log/WGS_GATK_$(date +"%Y%m%d%H%M").log"

set -ex
exec 2>>$logfile


# --------------------------
# paths to reference
# --------------------------
BUNDLE_RESOURCE=/N/slate/jiaji/CMG_references/hg38/GATK4_bundle_resource
REF=$BUNDLE_RESOURCE/Homo_sapiens_assembly38.fasta
annovar_db=/N/slate/jiaji/CMG_references/annovar/hg38/

DBSNP=$BUNDLE_RESOURCE/Homo_sapiens_assembly38.dbsnp138.vcf
KNOWNINDEL=$BUNDLE_RESOURCE/Homo_sapiens_assembly38.known_indels.vcf.gz
HAPMAP=$BUNDLE_RESOURCE/hapmap_3.3.hg38.vcf.gz
OMNI=$BUNDLE_RESOURCE/1000G_omni2.5.hg38.vcf.gz
G1000=$BUNDLE_RESOURCE/1000G_phase1.snps.high_confidence.hg38.vcf.gz
MILLS=$BUNDLE_RESOURCE/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
AXIOMPOLY=$BUNDLE_RESOURCE/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz


export OMP_NUM_THREADS=4

# --------------------------
# Functions
# --------------------------
function die() {
    echo "ERROR: $*" >&2
    exit 1
}

function echo_step() {
	echo "########################################"
    echo -e "# $(date '+%m-%d-%Y %H:%M:%S') - $1 "
    echo "########################################"
}

# --------------------------
# VQSR on SNP
# --------------------------

start_time_mod=$(date +%s)

## changed max-gaussians to 8 to match the Sentieon parameter
if [ ! -s $WD/WGS.GATK.joint.SNP.vqsr ]; then
	echo_step "Building SNP Recalibration Model..."
	gatk --java-options "-Xmx32g -Xms32g" \
		VariantRecalibrator \
		-R $REF \
		-V $JOINTVCF \
		--trust-all-polymorphic \
		-mode SNP \
		--max-gaussians 8 \
		-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
		--resource:hapmap,known=false,training=true,truth=true,prior=15 $HAPMAP \
		--resource:omni,known=false,training=true,truth=true,prior=12 $OMNI \
		--resource:1000G,known=false,training=true,truth=false,prior=10 $G1000 \
		--resource:dbsnp,known=true,training=false,truth=false,prior=7 $DBSNP \
		-O $WD/WGS.GATK.joint.SNP.vqsr \
		--tranches-file $WD/WGS.GATK.joint.SNP.tranches
fi

# apply the VQSR on SNP
if [ ! -s $WD/WGS.GATK.joint.SNP.vqsr.vcf.gz ]; then
	echo_step "Applying SNP Recalibration..."
	gatk --java-options "-Xmx32g -Xms32g" \
		ApplyVQSR \
		-V $JOINTVCF \
		-mode SNP \
		--recal-file $WD/WGS.GATK.joint.SNP.vqsr \
		--tranches-file $WD/WGS.GATK.joint.SNP.tranches \
		--truth-sensitivity-filter-level 99.7 \
		--create-output-variant-index true \
		-O $WD/WGS.GATK.joint.SNP.vqsr.vcf.gz 
fi

# --------------------------
# VQSR on INDEL
# --------------------------

## VQSR
## may need to reduce "--max-gaussians" (default 8 for SNP and 4 for Indel in Sentieon)
if [ ! -s $WD/WGS.GATK.joint.INDEL.vqsr ]; then
	echo_step "Building INDEL Recalibration Model..."
	gatk --java-options "-Xmx32g -Xms32g" \
		VariantRecalibrator \
		-V $WD/WGS.GATK.joint.SNP.vqsr.vcf.gz \
		-R $REF \
		--trust-all-polymorphic \
		-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
		-mode INDEL \
		--max-gaussians 4 \
		--resource:mills,known=false,training=true,truth=true,prior=12 $MILLS \
		--resource:axiomPoly,known=false,training=true,truth=false,prior=10 $AXIOMPOLY \
		--resource:dbsnp,known=true,training=false,truth=false,prior=2 $DBSNP \
		-O $WD/WGS.GATK.joint.INDEL.vqsr \
		--tranches-file $WD/WGS.GATK.joint.INDEL.tranches
fi

# apply the VQSR on INDEL
if [ ! -s $WD/WGS.GATK.joint.vqsr.vcf.gz ]; then
	echo_step "Building INDEL Recalibration Model..."
	gatk --java-options "-Xmx32g -Xms32g" \
		ApplyVQSR \
		-V $WD/WGS.GATK.joint.SNP.vqsr.vcf.gz \
		-mode INDEL \
		--recal-file $WD/WGS.GATK.joint.INDEL.vqsr \
		--tranches-file $WD/WGS.GATK.joint.INDEL.tranches \
		--truth-sensitivity-filter-level 99.7 \
		--create-output-variant-index true \
		-O $WD/WGS.GATK.joint.vqsr.vcf.gz 
fi


#filter the vcf file
if [ ! -s $WD/WGS.GATK.joint.vqsr.pass.vcf.gz ];then
    echo_step "Filter the vcf file..."
	vcftools --gzvcf $WD/WGS.GATK.joint.vqsr.vcf.gz --remove-filtered-all --recode --stdout | bgzip -@ 4 -c \
	> $WD/WGS.GATK.joint.vqsr.pass.vcf.gz || die "Filtering VCF failed"
fi

end_time_mod=$(date +%s)

if [ "$OSTYPE" = "darwin"* ]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
if [ "$OSTYPE" = "darwin"* ]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
echo "Module VQSR Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec">>$logfile


#start_time_mod=$(date +%s)

#variant annotation
#if [ ! -s $WD_anno/WGS.GATK.annovar.hg38_multianno.csv ];then
#    echo_step "Annotating the vcf file..."
#    convert2annovar.pl -format vcf4old \
#        $WD/WGS.GATK.joint.vqsr.pass.vcf.gz \
#        > $WD_anno/WGS.GATK.avinput || die "Conversion to annovar format failed"

#    table_annovar.pl \
#        $WD_anno/annovar/WGS.GATK.avinput \
#        $annovar_db \
#        -buildver hg38 \
#        -out $WD_anno/WGS.GATK.annovar \
#        -remove \
#        -protocol refGene,cytoBand,exac03,avsnp150,dbnsfp42c,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eur,gnomad312_genome,gnomad211_exome \
#        -operation g,r,f,f,f,f,f,f,f,f,f \
#        -nastring . -csvout -polish || die "Annotation failed"
#fi



#end_time_mod=$(date +%s)
#if [ "$OSTYPE" = "darwin"* ]; then start_date=$(date -j -f "%s" $start_time_mod); else start_date=$(date -d @$start_time_mod); fi
#if [ "$OSTYPE" = "darwin"* ]; then end_date=$(date -j -f "%s" $end_time_mod); else end_date=$(date -d @$end_time_mod); fi
#echo "Annotation Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec">>$logfile


# --------------------------
# HardFilter the SNP and INDEL (optional)
# --------------------------
# < 1hr
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
echo "Annotation Started: "$start_date"; Ended: "$end_date"; Elapsed time: "$(($end_time_mod - $start_time_mod))" sec">>$logfile

