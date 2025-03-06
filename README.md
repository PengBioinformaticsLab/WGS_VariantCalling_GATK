# WGS Variant Calling Pipeline

## 1. Merge BAM Files from Two Sequencing Runs (Optional)

```bash
sh . /scripts/00.all_sample_preprocess_merge.sh yourWorkDir/ sample.list
```
**Description:** Merges BAM files for samples from different sequencing runs.

**Estimated Runtime:** ~2 hours

**Output Directory:** `$WD/mergedBAMs/`

**Script:** `preprocess_merge_BAMs.sh`

**Command:**
```bash
java -Xmx24G -jar $picard MergeSamFiles \
    -I $WORKDIR/FirstRun/${SAMPLE}.hg38_bwa.sorted.bam \
    -I $WORKDIR/SecondRun/${SAMPLE}.hg38_bwa.sorted.bam \
    -O $WORKDIR/mergedBAMs/${SAMPLE}.sorted.bam --CREATE_INDEX true
```

---

## 2. Mark Duplicates and Variant Calling

```bash
sh ./scripts/01.all_sample_markdup_variantCall.sh yourWorkDir/ sample.list
```

### **MarkDuplicates**
**Estimated Runtime:** ~4 hours

**Script:** `markdup_bqsr_variantCall.sh`

**Command:**
```bash
gatk --java-options "-Xmx32g" \
    MarkDuplicatesSpark -I $WORKDIR/mergedBAMs/${SAMPLE}.sorted.bam \
    -O $WORKDIR/output/${SAMPLE}.sort.dup.bam \
    -M $WORKDIR/table/${SAMPLE}.marked_dup_metrics.txt \
    --create-output-bam-index true
```

### **Base Quality Score Recalibration (BQSR)**
**Estimated Runtime:** 6-10 hours

**BaseRecalibrator:**
```bash
gatk --java-options "-Xmx24g" BaseRecalibrator \
    -I $WORKDIR/output/${SAMPLE}.sort.dup.bam \
    -R $REF --known-sites $DBSNP --known-sites $KNOWNINDEL --known-sites $MILLS \
    -O $WORKDIR/table/${SAMPLE}.recalBefore.table
```

**ApplyBQSR:**
```bash
gatk ApplyBQSR \
    -R $REF \
    -I $WORKDIR/output/${SAMPLE}.sort.dup.bam \
    --bqsr-recal-file $WORKDIR/table/${SAMPLE}.recalBefore.table \
    -O $WORKDIR/output/bqsr/${SAMPLE}.sort.dup.bqsr.bam
```

### **HaplotypeCaller**
**Estimated Runtime:** ~24 hours

**Output Directory:** `variant_calls/gvcf/`

**Command:**
```bash
mkdir -p $WORKDIR/variant_calls/gvcf
gatk --java-options "-Xmx32g" HaplotypeCaller \
    -I $WORKDIR/output/bqsr/${SAMPLE}.sort.dup.bqsr.bam \
    -R $REF \
    -ERC GVCF \
    -O $WORKDIR/variant_calls/gvcf/${SAMPLE}.g.vcf.gz
```

---

## 3. Joint Genotyping

```bash
sh ./scripts/02.all_chr_jointgenotyping.sh yourWorkDir/ sample_gvcfs.list chromosome.list
```

### **GenomicsDBImport**
**Output Directory:** `$WORKDIR/GVCF_DB/${chr}_DB`

**Note:** Ensure WGS sample names match RNA-seq sample names; set column 1 as sample names in `sample_gvcfs.list`.

**Script:** `joint_genotyping.sh`

**Command:**
```bash
gatk --java-options "-Xmx32g -XX:+UseParallelGC -XX:ParallelGCThreads=10" \
    GenomicsDBImport \
    --genomicsdb-workspace-path $WORKDIR/GVCF_DB/${chr}_DB \
    --sample-name-map $WORKDIR/sample_gvcfs.list \
    -L $chr \
    --reader-threads 4
```

### **Joint Genotyping**
**Output Directory:** `$WORKDIR/variant_calls/vcf/${chr}.vcf.gz`

**Command:**
```bash
gatk --java-options "-Xmx32g" GenotypeGVCFs \
    -R $REF \
    -V gendb://$WORKDIR/GVCF_DB/${chr}_DB \
    -L $chr \
    -O $WORKDIR/variant_calls/vcf/${chr}.vcf.gz \
    --tmp-dir $WORKDIR/temp
```

---

## 4. Merge Chromosomes

**Output File:** `$WORKDIR/variant_calls/merged.vcf.gz`

**Script:** `merge_chrs.sh`

**Command:**
```bash
java -Xmx32G -jar $picard GatherVcfs \
    -I $WORKDIR/variant_calls/vcf/chr*.vcf.gz \
    -R $REF \
    -O $WORKDIR/variant_calls/merged.vcf.gz
```

---

## Notes
- Ensure `sample_gvcfs.list` contains the correct sample-to-file mapping.
- Adjust memory and thread settings based on available resources.

This completes the WGS variant calling pipeline using GATK. 