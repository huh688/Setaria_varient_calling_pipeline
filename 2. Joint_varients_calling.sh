#!/bin/bash
#SBATCH -p supermem
#SBATCH -t 7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%A.out

module load gatk/4.1.4.0

#reference genome file
REF=/scratch/haohu/index/bwa/Sviridis_500_v2.0

GenomicsDB=Sv_Diversity_446_DB_supermem

VCFlist=Sv_Diversity_446.list

TEMPDIR=/scratch/haohu/Diversity/temp_supermem

Intervals=/scratch/haohu/Diversity/Sviridis_500_v2.0.interval_list

Output=Sv_Diversity_446_JointCall_supermem

#Consolidate GVCFs with GenomicsDBImport
gatk GenomicsDBImport --genomicsdb-workspace-path $GenomicsDB -V $VCFlist --tmp-dir=$TEMPDIR -L $Intervals  --batch-size 5 --max-num-intervals-to-import-in-parallel 3

#Joint-Call Cohort with GenotypeGVCFs
gatk GenotypeGVCFs -R ${REF}.fa -V gendb://${GenomicsDB} -O ${Output}.vcf.gz --tmp-dir=$TEMPDIR

# Select SNPs and INDELs

gatk SelectVariants -R ${REF}.fa -V ${Output}.vcf.gz --select-type-to-include SNP -O ${Output}.SNP.vcf.gz

gatk SelectVariants -R ${REF}.fa -V ${Output}.vcf.gz --select-type-to-include INDEL -O ${Output}.INDEL.vcf.gz

# Hard filter for low quality varients - V2 - missing-values-evaluate-as-failing = TRUE
gatk VariantFiltration -R ${REF}.fa -V ${Output}.SNP.vcf.gz -O ${Output}.SNP.filtered.vcf.gz --filter-name "QD_FAIL" --filter-expression "QD < 2.00" --filter-name "SOR_FAIL" --filter-expression "SOR > 3.00" --filter-name "MQ_FAIL" --filter-expression "MQ < 40.00" --filter-name "MQRankSum_FAIL" --filter-expression "MQRankSum < -12.500" --filter-name "ReadPosRankSum_FAIL" --filter-expression "ReadPosRankSum < -8.000" --missing-values-evaluate-as-failing

gatk VariantFiltration -R ${REF}.fa -V ${Output}.INDEL.vcf.gz -O ${Output}.INDEL.filtered.vcf.gz --filter-name "QD_FAIL" --filter-expression "QD < 2.00" --filter-name "FS_FAIL" --filter-expression "FS > 200.000" --filter-name "ReadPosRankSum_FAIL" --filter-expression "ReadPosRankSum < -20.00" --filter-name "SOR_FAIL" --filter-expression "SOR > 10.00" --missing-values-evaluate-as-failing

# Select the varients passed filter - "PASS" only
gatk SelectVariants -R ${REF}.fa -V ${Output}.SNP.filtered.vcf.gz --exclude-filtered -O ${Output}.SNP.filtered.passed.vcf.gz

gatk SelectVariants -R ${REF}.fa -V ${Output}.INDEL.filtered.vcf.gz --exclude-filtered -O ${Output}.INDEL.filtered.passed.vcf.gz

# Merge SNP and INDEL vcfs 
gatk MergeVcfs -I ${Output}.SNP.filtered.passed.vcf.gz -I ${Output}.INDEL.filtered.passed.vcf.gz -O ${Output}.SNP.INDEL.filtered.passed.vcf.gz

gatk VariantFiltration -R ${REF}.fa -V Sv_Diversity.SNP.INDEL.filtered.passed.anno.vcf.gz -O Sv_Diversity.SNP.INDEL.filtered.passed.anno.photoperiod.vcf.gz -L Sv_Photoperiod.bed










