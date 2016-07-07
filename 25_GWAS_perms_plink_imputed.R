#################################################################
## GWAS permutations                                           ##
## by: Heather Wheeler 2016-01-28                              ##
## shuffle phenotypes and rerun GWAS n times                   ##
## output matrix of numSNPs x nPerms                           ##
#################################################################
args <- commandArgs(trailingOnly=T)
"%&%" <- function(a,b) paste(a,b,sep="")
date <- Sys.Date()
library(dplyr)

#directories
gwas.dir <- "/group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/"
dos.dir <- gwas.dir %&% "genotypes/UMich_imputation_results/mach_dosage_files/"
gt.dir <- gwas.dir %&% "genotypes/"

permset <- args[1]
#num perms
n <- 10
phenotype <- "rnGM412"
covariates <- "ageaudio,CumlCispdose,PC1-PC10"
covtag <- "age.cisp.10PCs"

#good SNP list
snpfile <- dos.dir %&% "N88.imputed_MAF0.05_INFOlt1.05_SNPlist.rds"
snplist <- as.character(readRDS(snpfile))

#results data.frame
resdf <- data.frame(SNP=snplist) %>% mutate(SNP=as.character(SNP))

#read in phenotypes and covariates
pheno <- read.table(gt.dir %&% "PtStudy.phenofile",header=T)
#pull FID,IID, and PCs, which will remain connected to genotypes
phenofixed <- dplyr::select(pheno,FID,IID,starts_with("PC"))
#pull subset of pheno with rnGM412 and covariates (set used in GWAS)
pheno2shuffle <- dplyr::filter(pheno,is.na(rnGM412)==FALSE,is.na(CumlCispdose)==FALSE,is.na(ageaudio)==FALSE) %>% dplyr::select(FID,rnGM412,ageaudio,CumlCispdose)

set.seed(as.numeric(permset)*123)

for(i in c(1:n)){
  #shuffle pheno & connect phenoperm to the wrong FID for merging with PCs
  phenoperm <- dplyr::sample_n(pheno2shuffle[,2:4],size=dim(pheno2shuffle)[1]) %>% mutate(FID=pheno2shuffle[,1])
  newpheno <- left_join(phenofixed,phenoperm,by="FID")
  write.table(newpheno, "tmp.pheno." %&% permset, quote=F, row.names=F)
  #run gwas
  startGWASout <- "cat " %&% gwas.dir %&% "GWAS_results/GWAS_dosage_results_header > tmp.gwas.permset" %&% permset
  system(startGWASout)
  for(j in c(1:22)){
    runECHO <- "echo " %&% dos.dir %&% "N88.imputed_maf0.01_R20.8_1000G.chr" %&% j %&% ".SNPxID.dose.gz " %&% dos.dir %&% "dose.list > plink-doselist.txt"
    runPLINK <- "plink --fam " %&% gt.dir %&% "N88.forImputation.fam --dosage plink-doselist.txt list sepheader format=1 --pheno tmp.pheno." %&% 
    permset %&% " --linear --pheno-name " %&% phenotype %&% " --covar tmp.pheno." %&% permset %&% " --covar-name " %&% covariates %&% " --map " %&% 
    dos.dir %&% "N88.imputed_maf0.01_R20.8_1000G.chr" %&% j %&% ".bim --out tmp.plink." %&% permset %&% "_chr" %&% j
    add2res1 <- "tail -n +2 tmp.plink." %&% permset %&% "_chr" %&% j %&% ".assoc.dosage > o." %&% permset
    add2res2 <- "cat o." %&% permset %&% " >> tmp.gwas.permset" %&% permset
    system(runECHO)
    system(runPLINK)
    system(add2res1)
    system(add2res2)
  }
  res <- read.table("tmp.gwas.permset" %&% permset,header=T) %>% mutate(SNP=as.character(SNP))
  #only include SNPs with maf>0.05 and info score < 1.05
  singleres <- dplyr::filter(res,SNP %in% snplist) %>% dplyr::select(P)
  colnames(singleres) <- c("P" %&% i)
  #add p-values to resdf
  resdf <- cbind(resdf, singleres)
}

write.table(resdf, "GWAS_permutations/N88.imputed_" %&% phenotype %&% "_" %&% covtag %&% "_" %&% n %&% "perms_set" %&% permset %&% ".txt",quote=F,row.names=F)
