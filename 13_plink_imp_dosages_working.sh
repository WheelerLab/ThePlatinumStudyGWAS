#!/bin/bash
#PBS -N plink
#PBS -S /bin/bash
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err

cd $PBS_O_WORKDIR
module load plink

##phenotype: rnGM412 covariates: ageaudio, new (3/28/16) CumlCispdose,PC1-10                                                                                                                                                                    
pheno="rnGM412"
covtag="age.NEWcisp.PCs1-10"
cat GWAS_results/GWAS_dosage_results_header > GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage
for (( i = 1 ; i <= 22; i++))
do
    echo /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose.gz /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/dose.list > plink-doselist-int.txt
    plink --fam genotypes/N88.forImputation.fam --dosage plink-doselist-int.txt list sepheader format=1 --pheno genotypes/PtStudy.phenofile_20160328 --linear --pheno-name ${pheno} --covar genotypes/PtStudy.phenofile_20160328 --covar-name ageaudio,CumlCispdose,PC1-PC10 --map genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim --out GWAS_results/N88.imputed_${pheno}_${covtag}_chr${i}
    tail -n +2 GWAS_results/N88.imputed_${pheno}_${covtag}_chr${i}.assoc.dosage > o
    cat o >> GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage
done
sort -gk10  GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage >  GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage.sorted

