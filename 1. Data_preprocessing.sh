#!/bin/bash
#SBATCH -p batch
#SBATCH -t 5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=hao.hu@okstate.edu
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%A.out

# GATK Pipeline 
# Version 2.1
# From FASTQ to BQSRed gVCF and BAM
# This is a  template, use "generate_individual_sbatch_from_gz.sh" for formatting sample-specific script, and then mutiple sample in parallel. 
# Authored by Hao Hu
# Last modifed: 09/23/2020
# Changes of V2.1 on 09/23/2020: Turn on all types of email notification; change output filename to "JobName.JobID.out"

#path to working directory
path=/scratch/haohu/Diversity

module load bwa/0.7.17
module load gatk/4.1.4.0
module load R/3.6.1

cd $path

# Type input file here
RawInput=TargetFileNameToChange

# Unzip GZIP file
/home/haohu/bin/unpigz -k ${RawInput}.fastq.gz

# Extract Read 1 and 2 from raw FASTQ file
python /home/haohu/bin/splitFASTQ.py ${RawInput}.fastq

read1=${RawInput}.fastq_R1.fastq

read2=${RawInput}.fastq_R2.fastq

# #reference genome file
REF=/scratch/haohu/index/bwa/Sviridis_500_v2.0

# BWA mapping to reference genome
bwa mem -M -t 32 $REF $read1 $read2 >${RawInput}.sam

# Add ReadGroups and convert SAM to BAM file
gatk AddOrReplaceReadGroups -I ${RawInput}.sam -O ${RawInput}.rg.sorted.bam -SO coordinate --RGID `head -1 $read1|awk -F: '{print $3":"$NF}'` --RGLB $RawInput --RGPU `head -1 $path/$read1|awk -F: '{print $3":"$NF}'` --RGPL illumina --RGSM $RawInput

# Mark duplicates in BAM
gatk MarkDuplicates -I ${RawInput}.rg.sorted.bam -O ${RawInput}.rg.sorted.dedup.bam --METRICS_FILE ${RawInput}.rg.sorted.dedup.bam.metrics

# Build BAM index
gatk BuildBamIndex -I ${RawInput}.rg.sorted.dedup.bam

# Run 1st round HaplotypeCaller to call raw VCF for BQSR
gatk HaplotypeCaller -R ${REF}.fa -I ${RawInput}.rg.sorted.dedup.bam -O ${RawInput}.rg.sorted.dedup.vcf.gz --native-pair-hmm-threads 32 

RawVCFInput=${RawInput}.rg.sorted.dedup

# Select SNPs and INDELs

gatk SelectVariants -R ${REF}.fa -V ${RawVCFInput}.vcf.gz --select-type-to-include SNP -O ${RawVCFInput}.SNP.vcf.gz

gatk SelectVariants -R ${REF}.fa -V ${RawVCFInput}.vcf.gz --select-type-to-include INDEL -O ${RawVCFInput}.INDEL.vcf.gz


# Hard filter for low quality varients
gatk VariantFiltration -R ${REF}.fa -V ${RawVCFInput}.SNP.vcf.gz -O ${RawVCFInput}.SNP.filtered.vcf.gz --filter-name "QD_FAIL" --filter-expression "QD < 2.00" --filter-name "SOR_FAIL" --filter-expression "SOR > 3.00" --filter-name "MQ_FAIL" --filter-expression "MQ < 40.00" --filter-name "MQRankSum_FAIL" --filter-expression "MQRankSum < -12.500" --filter-name "ReadPosRankSum_FAIL" --filter-expression "ReadPosRankSum < -8.000" --missing-values-evaluate-as-failing

gatk VariantFiltration -R ${REF}.fa -V ${RawVCFInput}.INDEL.vcf.gz -O ${RawVCFInput}.INDEL.filtered.vcf.gz --filter-name "QD_FAIL" --filter-expression "QD < 2.00" --filter-name "FS_FAIL" --filter-expression "FS > 200.000" --filter-name "ReadPosRankSum_FAIL" --filter-expression "ReadPosRankSum < -20.00" --filter-name "SOR_FAIL" --filter-expression "SOR > 10.00" --missing-values-evaluate-as-failing

# Select the varients passed filter - "PASS" only
gatk SelectVariants -R ${REF}.fa -V ${RawVCFInput}.SNP.filtered.vcf.gz --exclude-filtered -O ${RawVCFInput}.SNP.filtered.passed.vcf.gz

gatk SelectVariants -R ${REF}.fa -V ${RawVCFInput}.INDEL.filtered.vcf.gz --exclude-filtered -O ${RawVCFInput}.INDEL.filtered.passed.vcf.gz

# Merge SNP and INDEL vcfs 
gatk MergeVcfs -I ${RawVCFInput}.SNP.filtered.passed.vcf.gz -I ${RawVCFInput}.INDEL.filtered.passed.vcf.gz -O ${RawVCFInput}.SNP.INDEL.filtered.passed.vcf.gz

# Base Quality Score Recalibration, produce BQSR table from raw BAM
gatk BaseRecalibrator -R ${REF}.fa --known-sites ${RawVCFInput}.SNP.INDEL.filtered.passed.vcf.gz -I ${RawVCFInput}.bam -O ${RawVCFInput}.BQSR.table

# Apply BQSR table to raw BAM
gatk ApplyBQSR -R ${REF}.fa -I ${RawVCFInput}.bam --bqsr-recal-file ${RawVCFInput}.BQSR.table -O ${RawVCFInput}.BQSR.bam

# Get BQSR table from post-BQSR-BAM
gatk BaseRecalibrator -R ${REF}.fa --known-sites ${RawVCFInput}.SNP.filtered.passed.vcf.gz -I ${RawVCFInput}.BQSR.bam -O ${RawVCFInput}.postBQSR.table

# Get BQSR report (raw vs 1st round)

gatk AnalyzeCovariates -before ${RawVCFInput}.BQSR.table -after ${RawVCFInput}.postBQSR.table -plots ${RawVCFInput}.BQSRplot.pdf

# 2nd round HaplotypeCaller call, with BQSRed-BAM, produce GVCF
gatk HaplotypeCaller -R ${REF}.fa -I ${RawVCFInput}.BQSR.bam -O ${RawVCFInput}.BQSR.gvcf.gz --native-pair-hmm-threads 32 -ERC GVCF
