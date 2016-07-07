#!/bin/bash
#PBS -N R.gcta
#PBS -S /bin/bash
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err

cd $PBS_O_WORKDIR
module load R
module load gcta/1.24.4

gcta64 --bfile genotypes/gctaQC/N88_Recluster_TOP_20150911_FinalReport.postQC --autosome --make-grm-bin --maf 0.05 --out GCTA/N88_Recluster_TOP_20150911_FinalReport.postQC

gcta64 --reml --grm-bin GCTA/N88_Recluster_TOP_20150911_FinalReport.postQC --pheno GCTA/height_cm.phen  --out GCTA/N88_Recluster_TOP_20150911_FinalReport.postQC_height

gcta64 --reml --grm-bin GCTA/N88_Recluster_TOP_20150911_FinalReport.postQC --pheno GCTA/shuffled_height_cm.phen  --out GCTA/N88_Recluster_TOP_20150911_FinalReport.postQC_shuffled_height

gcta64 --reml --grm-bin GCTA/N88_Recluster_TOP_20150911_FinalReport.postQC --pheno GCTA/rnGM412.phen --out GCTA/N88_Recluster_TOP_20150911_FinalReport.postQC_rnGM412 --qcovar GCTA/ageaudio_CumlCispdose.qcovar

gcta64 --reml --grm-bin GCTA/N88_Recluster_TOP_20150911_FinalReport.postQC --pheno GCTA/ordCIPN8.phen --out GCTA/N88_Recluster_TOP_20150911_FinalReport.postQC_ordCIPN8 --qcovar GCTA/agediagnosis.qcovar
