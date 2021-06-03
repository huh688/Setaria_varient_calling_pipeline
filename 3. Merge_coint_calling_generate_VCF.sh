#!/bin/bash
#SBATCH -p batch
#SBATCH -t 5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%A.out

module load gatk/4.1.4.0

gatk MergeVcfs -I Sv_Diversity_446_JointCall_Chr1.vcf.gz  -I Sv_Diversity_446_JointCall_Chr2.vcf.gz -I Sv_Diversity_446_JointCall_Chr3.vcf.gz -I Sv_Diversity_446_JointCall_Chr4.vcf.gz -I Sv_Diversity_446_JointCall_Chr5.vcf.gz -I Sv_Diversity_446_JointCall_Chr6.vcf.gz -I Sv_Diversity_446_JointCall_Chr7.vcf.gz -I Sv_Diversity_446_JointCall_Chr8.vcf.gz -I Sv_Diversity_446_JointCall_Chr9.vcf.gz -I Sv_Diversity_446_JointCall_Chr10.vcf.gz -O Sv_Diversity_446_JointCall_merged.vcf.gz
