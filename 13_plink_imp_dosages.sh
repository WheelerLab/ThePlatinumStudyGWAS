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

##ototox phenotype: rnGM412 covariates: ageaudio, CumlCispdose                                                                                                                                                                                    
#covtag="age"
#cat GWAS_results/GWAS_dosage_results_header > GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage
#for (( i = 1 ; i <= 22; i++))
#do
#    echo /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose.gz /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/dose.list > plink-doselist.txt
#    plink --fam genotypes/N88.forImputation.fam --dosage plink-doselist.txt list sepheader format=1 --pheno genotypes/PtStudy.phenofile --linear --pheno-name rnGM412 --covar genotypes/PtStudy.covfile --covar-name ageaudio --map genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim --out GWAS_results/N88.imputed_rnGM412_${covtag}_chr${i}
#    tail -n +2 GWAS_results/N88.imputed_rnGM412_${covtag}_chr${i}.assoc.dosage > o
#    cat o >> GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage
#done
#sort -gk10  GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage >  GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage.sorted

##ototox phenotype: rnGM412 covariates: ageaudio, CumlCispdose
#covtag="age.cisp"
#cat GWAS_results/GWAS_dosage_results_header > GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage
#for (( i = 1 ; i <= 22; i++))
#do
#    echo /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose.gz /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/dose.list > plink-doselist.txt
#    plink --fam genotypes/N88.forImputation.fam --dosage plink-doselist.txt list sepheader format=1 --pheno genotypes/PtStudy.phenofile --linear --pheno-name rnGM412 --covar genotypes/PtStudy.covfile --covar-name ageaudio,CumlCispdose --map genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim --out GWAS_results/N88.imputed_rnGM412_${covtag}_chr${i}
#    tail -n +2 GWAS_results/N88.imputed_rnGM412_${covtag}_chr${i}.assoc.dosage > o
#    cat o >> GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage
#done
#sort -gk10  GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage >  GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage.sorted

##ototox phenotype: rnGM412 covariates: ageaudio, CumlCispdose, 10PCs
#covtag="age.cisp.10PCs"
#cat GWAS_results/GWAS_dosage_results_header > GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage
#for (( i = 1 ; i <= 22; i++))
#do
#    echo /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose.gz /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/dose.list > plink-doselist.txt
#    plink --fam genotypes/N88.forImputation.fam --dosage plink-doselist.txt list sepheader format=1 --pheno genotypes/PtStudy.phenofile --linear --pheno-name rnGM412 --covar genotypes/PtStudy.covfile --map genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim --out GWAS_results/N88.imputed_rnGM412_${covtag}_chr${i}
#    tail -n +2 GWAS_results/N88.imputed_rnGM412_${covtag}_chr${i}.assoc.dosage > o
#    cat o >> GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage
#done
#sort -gk10  GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage >  GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage.sorted

###NEXT: pseudo GWAS of PC1
#covtag=PC1
#cat GWAS_results/GWAS_dosage_results_header > GWAS_results/N88.imputed_${covtag}_chr1-22.assoc.dosage
#for (( i = 1 ; i <= 22; i++))
#do
#    echo /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose.gz /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/dose.list > plink-doselist.txt
#    plink --fam genotypes/N88.forImputation.fam --dosage plink-doselist.txt list sepheader format=1 --pheno genotypes/PtStudy.covfile --linear --pheno-name PC1 --map genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim --out GWAS_results/N88.imputed_${covtag}_chr${i}
#    tail -n +2 GWAS_results/N88.imputed_${covtag}_chr${i}.assoc.dosage > o
#    cat o >> GWAS_results/N88.imputed_${covtag}_chr1-22.assoc.dosage
#done
#sort -gk10  GWAS_results/N88.imputed_${covtag}_chr1-22.assoc.dosage >  GWAS_results/N88.imputed_${covtag}_chr1-22.assoc.dosage.sorted

##ototox phenotype: rnGM412 covariates: ageaudio, 10PCs                                                                                               
#covtag="age.10PCs"                                                                                                                                                                     
#cat GWAS_results/GWAS_dosage_results_header > GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage                                                                               
#for (( i = 1 ; i <= 22; i++))                                                                                                                                                               
#do                                                                                                                                                                                          
#    echo /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose.gz /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/dose.list > plink-doselist.txt                                                                                
#    plink --fam genotypes/N88.forImputation.fam --dosage plink-doselist.txt list sepheader format=1 --pheno genotypes/PtStudy.phenofile --linear --pheno-name rnGM412 --covar genotypes/PtStudy.covfile --covar-name ageaudio,PC1-PC10 --map genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim --out GWAS_results/N88.imputed_rnGM412_${covtag}_chr${i}                 
#    tail -n +2 GWAS_results/N88.imputed_rnGM412_${covtag}_chr${i}.assoc.dosage > o                                                                                                          
#    cat o >> GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage                                                                                                                
#done                                                                                                                                                                                        
#sort -gk10  GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage >  GWAS_results/N88.imputed_rnGM412_${covtag}_chr1-22.assoc.dosage.sorted 

##ototox phenotype: logGM412 covariates: ageaudio, CumlCispdose, 10PCs 
covtag="age.cisp.10PCs"
cat GWAS_results/GWAS_dosage_results_header > GWAS_results/N88.imputed_logGM412_${covtag}_chr1-22.assoc.dosage
for (( i = 1 ; i <= 22; i++))
do
    echo /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose.gz /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/dose.list > plink-doselist.txt
    plink --fam genotypes/N88.forImputation.fam --dosage plink-doselist.txt list sepheader format=1 --pheno genotypes/PtStudy.phenofile --linear --pheno-name logGM412 --covar genotypes/PtStudy.covfile --covar-name ageaudio,CumlCispdose,PC1-PC10 --map genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim --out GWAS_results/N88.imputed_logGM412_${covtag}_chr${i}
    tail -n +2 GWAS_results/N88.imputed_logGM412_${covtag}_chr${i}.assoc.dosage > o
    cat o >> GWAS_results/N88.imputed_logGM412_${covtag}_chr1-22.assoc.dosage
done
sort -gk10  GWAS_results/N88.imputed_logGM412_${covtag}_chr1-22.assoc.dosage >  GWAS_results/N88.imputed_logGM412_${covtag}_chr1-22.assoc.dosage.sorted

##phenotype: meanCIPN8 covariates: agediagnosis 
pheno="meanCIPN8"
covtag="agediagnosis"
cat GWAS_results/GWAS_dosage_results_header > GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage
for (( i = 1 ; i <= 22; i++))
do
    echo /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose.gz /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/dose.list > plink-doselist-int.txt
    plink --fam genotypes/N88.forImputation.fam --dosage plink-doselist-int.txt list sepheader format=1 --pheno genotypes/PtStudy.phenofile --linear --pheno-name ${pheno} --covar genotypes/PtStudy.phenofile --covar-name agediagnosis --map genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim --out GWAS_results/N88.imputed_${pheno}_${covtag}_chr${i} 
    tail -n +2 GWAS_results/N88.imputed_${pheno}_${covtag}_chr${i}.assoc.dosage > o
    cat o >> GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage
done
sort -gk10  GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage >  GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage.sorted

##phenotype: dosegroup covariates: 10 PCs
pheno="dosegroup"
covtag="10PCs"
cat GWAS_results/GWAS_dosage_results_header > GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage
for (( i = 1 ; i <= 22; i++))
do
    echo /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose.gz /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/dose.list > plink-doselist-int.txt
    plink --fam genotypes/N88.forImputation.fam --dosage plink-doselist-int.txt list sepheader format=1 --pheno genotypes/PtStudy.phenofile --linear --pheno-name ${pheno} --covar genotypes/PtStudy.phenofile --covar-name PC1-PC10 --map genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim --out GWAS_results/N88.imputed_${pheno}_${covtag}_chr${i} --1
    tail -n +2 GWAS_results/N88.imputed_${pheno}_${covtag}_chr${i}.assoc.dosage > o
    cat o >> GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage
done
sort -gk10  GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage >  GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage.sorted

##phenotype: rnGM48norway covariates: ageaudio
pheno="rnGM48nor"
covtag="ageaudio"
phenfile="Norway.phenofile"

cat GWAS_results/GWAS_dosage_results_header > GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage
for (( i = 1 ; i <= 22; i++))
do
    echo /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose.gz /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/dose.list > plink-doselist-int.txt
    plink --fam genotypes/N88.forImputation.fam --dosage plink-doselist-int.txt list sepheader format=1 --pheno genotypes/${phenfile} --linear --pheno-name ${pheno} --covar genotypes/${phenfile} --covar-name ${covtag} --map genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim --out GWAS_results/N88.imputed_${pheno}_${covtag}_chr${i}
    tail -n +2 GWAS_results/N88.imputed_${pheno}_${covtag}_chr${i}.assoc.dosage > o
    cat o >> GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage
done
sort -gk10  GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage >  GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage.sorted

##ototox phenotype: GM412 covariates: ageaudio, CumlCispdose, 10PCs                                                                                                                                                
covtag="age.cisp.10PCs"
phentag="GM412"
cat GWAS_results/GWAS_dosage_results_header > GWAS_results/N88.imputed_${phentag}_${covtag}_chr1-22.assoc.dosage
for (( i = 1 ; i <= 22; i++))
do
    echo /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose.gz /group/dolan-lab/hwheeler/ThePlatinumStudy/GW\
AS/genotypes/UMich_imputation_results/mach_dosage_files/dose.list > plink-doselist.txt
    plink --fam genotypes/N88.forImputation.fam --dosage plink-doselist.txt list sepheader format=1 --pheno genotypes/PtStudy.phenofile --linear --pheno-name ${phentag} --covar genotypes/PtStudy.covfile --covar\
-name ageaudio,CumlCispdose,PC1-PC10 --map genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim --out GWAS_results/N88.imputed_${phentag}_${covtag}_chr${i}
    tail -n +2 GWAS_results/N88.imputed_${phentag}_${covtag}_chr${i}.assoc.dosage > o
    cat o >> GWAS_results/N88.imputed_${phentag}_${covtag}_chr1-22.assoc.dosage
done
sort -gk10  GWAS_results/N88.imputed_${phentag}_${covtag}_chr1-22.assoc.dosage >  GWAS_results/N88.imputed_${phentag}_${covtag}_chr1-22.assoc.dosage.sorted

##ototox phenotype: logGM412 covariates: ageaudio, CumlCispdose, 10PCs                                                                                                                                             
covtag="age.cisp.10PCs"
phentag="logGM412"
cat GWAS_results/GWAS_dosage_results_header > GWAS_results/N88.imputed_${phentag}_${covtag}_chr1-22.assoc.dosage
for (( i = 1 ; i <= 22; i++))
do
    echo /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose.gz /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/dose.list > plink-doselist.txt
    plink --fam genotypes/N88.forImputation.fam --dosage plink-doselist.txt list sepheader format=1 --pheno genotypes/PtStudy.phenofile --linear --pheno-name ${phentag} --covar genotypes/PtStudy.covfile --covar-name ageaudio,CumlCispdose,PC1-PC10 --map genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim --out GWAS_results/N88.imputed_${phentag}_${covtag}_chr${i}
    tail -n +2 GWAS_results/N88.imputed_${phentag}_${covtag}_chr${i}.assoc.dosage > o
    cat o >> GWAS_results/N88.imputed_${phentag}_${covtag}_chr1-22.assoc.dosage
done
sort -gk10  GWAS_results/N88.imputed_${phentag}_${covtag}_chr1-22.assoc.dosage >  GWAS_results/N88.imputed_${phentag}_${covtag}_chr1-22.assoc.dosage.sorted

##phenotype: meanCIPN8 covariates: agediagnosis                                                                                                                                                                    
pheno="meanCIPN8"
covtag="agediagnosis"
cat GWAS_results/GWAS_dosage_results_header > GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage
for (( i = 1 ; i <= 22; i++))
do
    echo /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose.gz /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/mach_dosage_files/dose.list > plink-doselist-int.txt
    plink --fam genotypes/N88.forImputation.fam --dosage plink-doselist-int.txt list sepheader format=1 --pheno genotypes/PtStudy.phenofile --linear --pheno-name ${pheno} --covar genotypes/PtStudy.phenofile --covar-name agediagnosis --map genotypes/UMich_imputation_results/mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim --out GWAS_results/N88.imputed_${pheno}_${covtag}_chr${i}
    tail -n +2 GWAS_results/N88.imputed_${pheno}_${covtag}_chr${i}.assoc.dosage > o
    cat o >> GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage
done
sort -gk10  GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage >  GWAS_results/N88.imputed_${pheno}_${covtag}_chr1-22.assoc.dosage.sorted

