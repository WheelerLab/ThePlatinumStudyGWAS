#!/bin/bash                                                                                                                                                               
#PBS -N vcf.4.impute
#PBS -S /bin/bash
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err

module load plink
module load vcftools
module load tabix/0.2.6

cd /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results

for ( i = 1 ; i <= 21; i++)
do
    ## use vcftools to remove indels and remove SNPs with MAF<0.01, keep INFO 
    #vcftools --gzvcf chr$i.dose.vcf.gz --remove-indels --maf 0.01 --recode --recode-INFO-all --stdout | gzip -c > chr$i.dose_maf0.01_rm.indel.vcf.gz
    #pull data of interest
    perl pull_qual_info.pl chr$i.dose_maf0.01_rm.indel.vcf.gz > chr$i.r2
done
##plot R2 "Estimated Imputation Accuracy" and ER2 (concordance) "Empirical (Leave-One-Out) R-square (available only for genotyped variants)"
R --vanilla < plot_impute_R2_ER2.r

##plink command with dose.gz files made in 12_vcf2mach.dosage_1000G.pl
plink --fam ../../N88.forImputation.fam --dosage c.txt list sepheader format=1 --pheno ../../PtStudy.phenofile --linear --pheno-name rnGM412 --out test.dosage --covar ../../PtStudy.covfile