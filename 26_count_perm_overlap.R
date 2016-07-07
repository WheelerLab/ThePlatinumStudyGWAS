#####################################################################
## Count overlap of GWAS permutations with LCL cytotox             ##
## by: Heather Wheeler 2016-02-01                                  ##
## count #SNPs less than some p-value threshold (pGWAS) that are   ##
## also associated with cytotox and some p-value threshold (pLCL)  ##
#####################################################################
args <- commandArgs(trailingOnly=T)
"%&%" <- function(a,b) paste(a,b,sep="")
date <- Sys.Date()
library(dplyr)

#directories
gwas.dir <- "/group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/"
perm.dir <- gwas.dir %&% "GWAS_permutations/" 
lcl.dir <- gwas.dir %&% "LCL_cytotoxicity/output/"
obs.dir <- gwas.dir %&% "GWAS_results/"
enrich.dir <- gwas.dir %&% "LCL_cytotoxicity/enrichment_results/"

#run once to make LCL cisplatin cytotoxicity GWAS results RDS file
#phenolist <- scan("/group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/LCL_cytotoxicity/pheno.list","c")
#for(pheno in phenolist){
#  lclfile <- lcl.dir %&% "CEU.TGP_and_imputed.rmBAD.20130526." %&% pheno %&% ".assoc.txt"
#  lcl <- read.table(lclfile,header=T)
#  shortlcl <- dplyr::select(lcl,rs,p_wald)
#  colnames(shortlcl) <- c("SNP", "P")
#  saveRDS(shortlcl,lcl.dir %&% "p_wald.CEU.TGP_and_imputed.rmBAD.20130526." %&% pheno %&% "assoc.RDS")
#}

l <- as.numeric(args[1]) #LCL p thresh
g <- as.numeric(args[2]) #GWAS p thresh
lclpheno <- args[3] #LCL drug phenotype
  
lcl <- readRDS(lcl.dir %&% "p_wald.CEU.TGP_and_imputed.rmBAD.20130526." %&% lclpheno %&% ".phenoassoc.RDS") #for faster reading

#topObs <- read.table(obs.dir %&% "N88.imputed_rnGM412_age.cisp.10PCs_chr1-22.assoc.dosage.sorted.P_0.05",header=T) #for imputed
topObs <- read.table(obs.dir %&% "N88_Recluster_TOP_20150911_FinalReport.forPCA_rnGM412_ageaudio.cisp.10PCs.assoc.linear.adjusted",header=T) %>% 
  mutate(P=UNADJ) #for genotyped


firstperm <- read.table(perm.dir %&% "N88.imputed_rnGM412_age.cisp.10PCs_10perms_set1.txt",header=T)
allperms <- dplyr::select(firstperm,SNP)
#read in all perm files and add to data.frame
for(i in c(1:100)){
  perm <- read.table(perm.dir %&% "N88.imputed_rnGM412_age.cisp.10PCs_10perms_set" %&% i %&% ".txt",header=T)
  newperm <- perm[unique(perm$SNP),]
  allperms <- cbind(allperms, perm[,2:11])
}
  
#function that calculates the empirical P value by calculating the overlap between LCL and each GWAS permutation and comparing to the observed overlap
#requires LCL P-value threshold and GWAS P-value threshold, default is not to plot a histogram
countOverlap <- function(lclP,gwasP,lclpheno,hist=FALSE,histfile=FALSE){
  #get LCL SNP list with P < lclP threshold
  lclList <- dplyr::filter(lcl,P < lclP) %>% dplyr::select(SNP) %>% mutate(SNP = as.character(SNP))
  #get GWAS SNP list with P < gwasP threshold
  obsList <- dplyr::filter(topObs, P < gwasP) %>% dplyr::select(SNP) %>% mutate(SNP = as.character(SNP))
  #calculate the observed number of overlap SNPs
  if(dim(table(obsList$SNP %in% lclList$SNP))==1){
    res <- c(lclpheno, lclP, gwasP, 0, NA)
    return(res)
  }else{
    obsOverlap <- table(obsList$SNP %in% lclList$SNP)[[2]]
    #count # times permuted overlap is greater than observed to get empirical P-value
    permcounts <- vector() #vector to hold perm counts
    for(i in c(2:dim(allperms)[2])){
      l <- allperms[allperms[,i] < gwasP,] #truncate all perms table to those SNPs w/P < gwasP in column i
      permList <- l$SNP #get perm i SNP list with P < gwasP threshold
      if(dim(table(permList %in% lclList$SNP))==1){
        permOverlap = 0
      }else{
        permOverlap <- table(permList %in% lclList$SNP)[[2]] #caculate number of LCL/perm i overlap SNPs
      }
      permcounts <- c(permcounts, permOverlap) #add count to permcounts vector
    }
    if(max(permcounts) > obsOverlap){ #if at least 1 permutation had more overlap than observed, calc empP
      empP <- table(permcounts > obsOverlap)[[2]]/(dim(allperms)[2]-1)
    }else{ #otherwise empP < 1/(# perms)
      empP <- 0 
    }  
    if(hist==TRUE){
      options(bitmapType='cairo') #for png to work on linux with no X11
      png(filename=histfile)
      hist(permcounts,xlab="Number of Overlap SNPs",main="LCL P<" %&% lclP %&% " & GWAS P<" %&% gwasP %&% " empP = " %&% empP, n=20)
      points(obsOverlap,0,pch=19,cex=3)
      dev.off()
    }
    res <- c(lclpheno, lclP, gwasP, obsOverlap, empP)
    return(res)
  }
}

permP <- countOverlap(l, g, lclpheno, hist=TRUE,histfile=enrich.dir %&% "GWASp."  %&%
  g %&% "_LCLp" %&% l %&% "_" %&% lclpheno %&% "_" %&% date %&% ".png")
#names(permP) <- c('lclpheno','lclP','gwasP','obsOverlap','empP')

#print empirical P-values
outfile <- enrich.dir %&% "PtStudy_gt_rnGM412_ageaudio.cisp.10PCs.assoc.linear_GWASp."  %&% 
  g %&% "_LCLp" %&% l %&% "_" %&% lclpheno %&% "_" %&% date %&% ".txt"
write(permP, outfile, ncolumns = 5)
