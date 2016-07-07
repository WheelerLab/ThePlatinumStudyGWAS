#!/bin/bash
#PBS -N gemma
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err

cd $PBS_O_WORKDIR

##Run gemma on imputed CEU CIS_log2IC50 phenotype
module load gemma/0.94

dir=/group/dolan-lab/PREVIOUS\ MEMBERS/WHEELER\,\ Heather\ 2010-2015/1KGP_imputation/

for pt in `cat LCL_cytotoxicity/pheno.list`;
do
    gemma -g "${dir}"bimbam_files/CEU.TGP_and_imputed.rmBAD.20130526.geno -p LCL_cytotoxicity/${pt} -gk 2 -o GRM.CEU.TGP_and_imputed.rmBAD.20130526.${pt}
    gemma -g "${dir}"bimbam_files/CEU.TGP_and_imputed.rmBAD.20130526.geno -p LCL_cytotoxicity/${pt} -k output/GRM.CEU.TGP_and_imputed.rmBAD.20130526.${pt}.sXX.txt -a "${dir}"bimbam_files/CEU.TGP_and_imputed.rmBAD.20130526.snp.info -lmm 4 -o CEU.TGP_and_imputed.rmBAD.20130526.${pt}
done
