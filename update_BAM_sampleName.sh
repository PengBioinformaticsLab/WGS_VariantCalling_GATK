#!/bin/bash
#SBATCH -A r00346
#SBATCH --mail-user=jiaji@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=24gb
#SBATCH --partition=general
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -J WGS_reheadSM
#SBATCH -o log/WGS_reheadSM_%j.out

# Usage: ./update_BAM_sampleName.sh sample_list.txt
OUTDIR=$1

while IFS=$'\t' read -r new_sm bam; do
    # Skip empty lines or headers
    [[ -z "$bam" || -z "$new_sm" ]] && continue

    echo "Processing $bam -> SM=$new_sm"

    # 1. Extract the BAM header
    header=$(mktemp)
    samtools view -H "$bam" > "$header"

    # 2. Update all @RG SM tags in the header
    new_header=$(mktemp)
    awk -v new_sm="$new_sm" '
        BEGIN {OFS="\t"}
        /^@RG/ {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^SM:/) {
                    $i = "SM:" new_sm
                }
            }
        }
        {print}
    ' "$header" > "$new_header"

    # 3. Reheader the BAM with the new SM tag
    output_bam=$OUTDIR/${new_sm}.hg38_bwa.sorted.bam
    samtools reheader -P "$new_header" "$bam" > "$output_bam"

    # 4. Index the new BAM
    samtools index "$output_bam"

    # Cleanup temporary files
    rm "$header" "$new_header"

done < "$1"
