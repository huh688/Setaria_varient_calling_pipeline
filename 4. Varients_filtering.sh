#!/bin/bash
#SBATCH -p batch
#SBATCH -t 5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%A.out

module load gatk/4.1.4.0

# Type input file here
RawInput=Sv_Diversity_446_JointCall_merged

# #reference genome file
REF=/scratch/haohu/index/bwa/Sviridis_500_v2.0

# Select SNPs and INDELs, disable --missing-values-evaluate-as-failing, enable --genotype-filter-expression "DP < 2 || DP > 50" --genotype-filter-name "DP_Fail"

gatk SelectVariants -R ${REF}.fa -V ${RawInput}.vcf.gz --select-type-to-include SNP -O ${RawInput}.SNP.vcf.gz

gatk VariantFiltration -R ${REF}.fa -V ${RawInput}.SNP.vcf.gz -O ${RawInput}.SNP.filtered.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" --genotype-filter-expression "DP < 2 || DP > 50" --genotype-filter-name "DP_Fail"

# Enable -restrict-alleles-to BIALLELIC
gatk SelectVariants -R ${REF}.fa -V ${RawInput}.SNP.filtered.vcf.gz --exclude-filtered -restrict-alleles-to BIALLELIC -O ${RawInput}.SNP.filtered.passed.vcf.gz

# Select INDELs, disable --missing-values-evaluate-as-failing

gatk SelectVariants -R ${REF}.fa -V ${RawInput}.vcf.gz --select-type-to-include INDEL -O ${RawInput}.INDEL.vcf.gz

gatk VariantFiltration -R ${REF}.fa -V ${RawInput}.INDEL.vcf.gz -O ${RawInput}.INDEL.filtered.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

gatk SelectVariants -R ${REF}.fa -V ${RawInput}.INDEL.filtered.vcf.gz --exclude-filtered -O ${RawInput}.INDEL.filtered.passed.vcf.gz

#Filter with vcftools

/home/haohu/bin/vcftools-0.1.16/bin/vcftools --gzvcf ${RawInput}.SNP.filtered.passed.vcf.gz --max-missing 0.9 --maf 0.01 --recode --recode-INFO-all --stdout | bgzip -c > ${RawInput}.SNP.filtered.passed.vcftools.MAXMISS09.MAF001.gz

tabix ${RawInput}.SNP.filtered.passed.vcftools.MAXMISS09.MAF001.gz