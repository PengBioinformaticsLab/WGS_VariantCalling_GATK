WORKDIR=$1
SAMPLE=$2


cd $WORKDIR
mkdir -p $WORKDIR/log/jobs

cat <<EOF > $WORKDIR/log/jobs/${SAMPLE}.mkdp.vc.job
#!/bin/bash

#SBATCH -A r00302
#SBATCH --mail-user=jiaji@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mem=64gb
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH -J ${SAMPLE}.mkdp.vc
#SBATCH -o log/${SAMPLE}_mkdp_vc_%j.out

export OMP_NUM_THREADS=8
sh $WORKDIR/scripts/markdup_bqsr_variantCall.sh $WORKDIR $SAMPLE

EOF

sbatch $WORKDIR/log/jobs/${SAMPLE}.mkdp.vc.job
