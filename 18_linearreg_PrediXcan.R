########################################
#### linear regression PrediXcan   ####
#### by Heather E. Wheeler 20151226 ####
########################################
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)

tis <- args[1]
phen <- args[2]

#for testing:
#tis <- "TW_Muscle-Skeletal_ElasticNet.0.5"
#phen <- "rnGM412"

#my.dir <- "/group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/"
my.dir <- "/group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/"
exp.dir <- my.dir %&% "genotypes/UMich_imputation_results/PrediXcan/pred_exp/"
pheno.dir <- my.dir %&% "genotypes/"
res.dir <- my.dir %&% "genotypes/UMich_imputation_results/PrediXcan/PrediXcan_results/"

expfile <- exp.dir %&% "N88_n953_" %&% tis %&% ".txt"
phenofile <- pheno.dir %&% "PtStudy.phenofile"
covfile <- pheno.dir %&% "PtStudy.covfile"
gencodefile <- exp.dir %&% "gencode.v18.genes.patched_contigs.summary.protein"

########functions
linreg <- function(gt,pt=phenotype,covariates=covmat){
  m<-summary(lm(pt ~ gt + covariates))
  n<-length(pt)-length(m$na.action)
  res <- c(m$coef[2,],n)
  names(res) <- c('Estimate','SE','t-value','P','n')
  return(res)
}

exp <- read.table(expfile,header=TRUE)
#transpose to be like genotypes
texp <- t(exp)
#for testing:
#texp <- texp[1:10,]
#only include genes with mean exp > 0
texp<-texp[rowMeans(texp)>0,]
dim(texp)


pheno <- read.table(phenofile,header=T)
cov <- read.table(covfile,header=T) %>% dplyr::select(FID,starts_with("PC"))
pheno <- left_join(pheno,cov,by="FID")

phenotype <- dplyr::select(pheno,starts_with(phen))[,1]

gencode <- read.table(gencodefile)
colnames(gencode) <- c("chr","strand","start","end","ensid","gene","type","status")
##run linear regression:

###1. covariates: ageaudio, CumlCispdose
covmat <- dplyr::select(pheno,ageaudio,CumlCispdose)
covmat <- as.matrix(covmat)
lingwas <- apply(texp,1,linreg)
if(tis == "DGN-WB"){
  tlingwas <- data.frame(t(lingwas)) %>% dplyr::mutate(gene=colnames(lingwas)) 
  output <- left_join(tlingwas,gencode,by="gene") %>% dplyr::select(-type,-status) %>% arrange(P)
  }else{
  tlingwas <- data.frame(t(lingwas)) %>% dplyr::mutate(ensid=colnames(lingwas)) 
  output <- left_join(tlingwas,gencode,by="ensid") %>% dplyr::select(-type,-status) %>% arrange(P)
}

write.table(output, file=res.dir %&% "N88_n953_" %&% phen %&% "_ageaudio.cisp_" %&% tis %&% ".linreg.PrediXcan"
            ,quote=F,row.names=F)

##2. covariates: ageaudio, CumlCispdose, PCs1-10
covmat <- dplyr::select(pheno,ageaudio,CumlCispdose,starts_with("PC"))
covmat <- as.matrix(covmat)
lingwas <- apply(texp,1,linreg)
if(tis == "DGN-WB"){
  tlingwas <- data.frame(t(lingwas)) %>% dplyr::mutate(gene=colnames(lingwas)) 
  output <- left_join(tlingwas,gencode,by="gene") %>% dplyr::select(-type,-status) %>% arrange(P)
}else{
  tlingwas <- data.frame(t(lingwas)) %>% dplyr::mutate(ensid=colnames(lingwas)) 
  output <- left_join(tlingwas,gencode,by="ensid") %>% dplyr::select(-type,-status) %>% arrange(P)
}

write.table(output, file=res.dir %&% "N88_n953_" %&% phen %&% "_ageaudio.cisp.10PCs_" %&% tis %&% ".linreg.PrediXcan"
            ,quote=F,row.names=F)
