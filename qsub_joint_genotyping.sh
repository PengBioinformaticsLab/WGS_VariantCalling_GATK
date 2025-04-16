WORKDIR=$1
sample_gvcf_map=$2
chr=$3

cd $WORKDIR
mkdir -p $WORKDIR/temp $WORKDIR/GVCF_DB $WORKDIR/variant_calls/vcf
rm -rf $WORKDIR/temp/*


cat <<EOF > $WORKDIR/log/jobs/${chr}.joint_genotyping.job
#!/bin/bash

#SBATCH -A r00302
#SBATCH --mail-user=jiaji@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=90:00:00
#SBATCH --mem=12gb
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH -J ${chr}_genotyping
#SBATCH -o log/${chr}_genotyping_%j.out


sh $WORKDIR/scripts/joint_genotyping.sh $WORKDIR $sample_gvcf_map $chr

EOF

sbatch $WORKDIR/log/jobs/${chr}.joint_genotyping.job
