#!/bin/bash                                                                                                                                                               
#PBS -N vcf.4.impute
#PBS -S /bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err

###################################################################                                                                                                      
#    THE PLATINUM STUDY RIKEN GWAS IMPUTATION WORKFLOW            #                                                                                                       
#    Heather E. Wheeler 2015-11-19                                #                                                                                                       
###################################################################   

module load plink
module load vcftools
module load tabix/0.2.6

cd /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/

for (( i = 1 ; i <= 23; i++)) 
do
    plink --bfile N88.forImputation --chr $i --recode vcf --out N88.forImputation.chr$i
    vcf-sort N88.forImputation.chr$i.vcf | bgzip -c > N88.forImputation.chr$i.vcf.gz
done

#chr23 is non-pseudoautosomal region of X, all male
plink --bfile N88.forImputation --chr 23 --recode vcf --out N88.forImputation.chr23
vcf-sort N88.forImputation.chr23.vcf | bgzip -c > N88.forImputation.chr23.vcf.gz
#imputation fail: No valid chromosomes found!, try replacing 23 with X
perl change23toX.pl N88.forImputation.chr23.vcf.gz  | bgzip -c > N88.forImputation.chrX.vcf.gz
#uploaded to UMich: "Chromosome X imputation is currently under preperation. Please upload only autosomale chromosomes." chrX not yet available on server


##use sftp to upload *vcf.gz files to https://imputationserver.sph.umich.edu/start.html#!pages/run, see 'vcf.filelist.for.sftp' for paths
##above didn't work, mounted tarbell to Desktop

##started QC run for all chromosomes, expecting strand flip errors
##Options selected in GUI: 
## Reference Panel:1000G 1 v3 Shapeit2 (no singletons)
## Phasing: SHAPEIT
## Population: EUR
## Mode: Quality Control Only

##download statistics.txt file from Imputationserver and renamed 'statistics.N88.chr1-22.txt' to get errors
##To fix strand flips:

#get snp list to flip strand
grep FILTER\ -\ Strand statistics.N88.chr1-22.txt > test
grep -o rs[0-9]* test > test1
grep -o exm[0-9]* test > test2
cat test1 test2 > N88.chr1-22.strand.switches

#get snp list to remove one duplicate
grep Duplicate statistics.N88.chr1-22.txt > test
grep -o rs[0-9]* test >test1
grep -o exm[0-9]* test > test2
cat test1 test2 > N88.chr1-22.dup.pos

for (( i = 1 ; i <= 22; i++))
do
    plink --bfile N88.forImputation --chr $i --exclude N88.chr1-22.dup.pos --flip N88.chr1-22.strand.switches --recode vcf --out N88.forImputation.chr$i.strand.switches
    vcf-sort N88.forImputation.chr$i.strand.switches.vcf | bgzip -c > N88.forImputation.chr$i.strand.switches.vcf.gz
done

##mount tarbell to upload *strand.switches.vcf.gz files to https://imputationserver.sph.umich.edu/start.html#!pages/run
## Use same options as above expect Mode: Quality Control and Imputation

#job failed, remove additional problem snps
##download statistics.txt file from Imputationserver and renamed 'statistics2.N88.chr1-22.txt' to get errors 
grep FILTER statistics2.N88.chr1-22.txt > test
grep -o rs[0-9]* test > test1
grep -o exm-rs[0-9]* test > test2
grep -o VG[0-9]*S[0-9]* test > test3
grep D genotypes/N88.forImputation.bim | cut -f 2 > test4 #rm indels
cat N88.chr1-22.dup.pos test1 test2 test3 test4 > N88.chr1-22.exclude.SNPlist

for (( i = 1 ; i <= 22; i++))
do
    plink --bfile N88.forImputation --chr $i --exclude N88.chr1-22.exclude.SNPlist --flip N88.chr1-22.strand.switches --recode vcf --out N88.forImputation.chr$i.strand.switches
    vcf-sort N88.forImputation.chr$i.strand.switches.vcf | bgzip -c > N88.forImputation.chr$i.strand.switches.vcf.gz
done



##save output in /group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/

##to get results:
##click each arrow to get private link on results page
wget https://imputationserver.sph.umich.edu/share/hwheeler/7be372baa8b8f848acecc4943d67339f/chr_5.zip
##etc. different link for each chr, couldn't figure out how to loop through all
unzip chr_11.zip #then type in password, couldn't figure out how to loop through all
#Archive:  chr_11.zip
#[chr_11.zip] chr11.dose.vcf.gz password: ##was emailed: MMhulrpPeP 


##filter results and check accuracy plots
for (( i = 1 ; i <= 22; i++))
do
    ## use vcftools to remove indels and remove SNPs with MAF<0.01, keep INFO 
    vcftools --gzvcf chr$i.dose.vcf.gz --remove-indels --maf 0.01 --recode --recode-INFO-all --stdout | gzip -c > chr$i.dose_maf0.01_rm.indel.vcf.gz
    #pull data of interest
    perl pull_qual_info.pl chr$i.dose_maf0.01_rm.indel.recode.vcf.gz > chr$i.r2
done
##plot R2 "Estimated Imputation Accuracy" and ER2 (concordance) "Empirical (Leave-One-Out) R-square (available only for genotyped variants)"
R --vanilla < plot_impute_R2_ER2.r
