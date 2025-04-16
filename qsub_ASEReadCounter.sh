WORKDIR=$1
SAMPLE=$2
BAM=$3
VCF=$4

cd $WORKDIR

cat <<EOF > $WORKDIR/log/jobs/${SAMPLE}.ase.job
#!/bin/bash

#SBATCH -A r00302
#SBATCH --mail-user=jiaji@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=32gb
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH -J ${SAMPLE}.ase.readcounter
#SBATCH -o log/${SAMPLE}.ase.%j.out

export OMP_NUM_THREADS=4
sh $WORKDIR/scripts/ASEreadcounter.sh $WORKDIR $SAMPLE $BAM $VCF

EOF

sbatch $WORKDIR/log/jobs/${SAMPLE}.ase.job

