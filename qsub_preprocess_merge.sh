WORKDIR=$1
SAMPLE=$2

BAM1=$WORKDIR/BAM/20240410_LH00300_0040_B22555FLT4/${SAMPLE}.hg38_bwa.sorted.bam
BAM2=$WORKDIR/BAM/20240503_LH00300_0043_A227KTFLT4/${SAMPLE}.hg38_bwa.sorted.bam

cd $WORKDIR
mkdir -p $WORKDIR/mergedBAMs
mkdir -p $WORKDIR/log/jobs

cat <<EOF > $WORKDIR/log/jobs/${SAMPLE}.preprocess_merge.job
#!/bin/bash

#SBATCH -A r00346
#SBATCH --mail-user=jiaji@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --mem=32gb
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH -J ${SAMPLE}_merge
#SBATCH -o log/${SAMPLE}_merge_%j.out


sh $WORKDIR/scripts/preprocess_merge_BAMs.sh $BAM1 $BAM2 $SAMPLE $WORKDIR

EOF

sbatch $WORKDIR/log/jobs/${SAMPLE}.preprocess_merge.job
