########################################
#### Ordinal regression PrediXcan   ####
#### by Heather E. Wheeler 20151226 ####
########################################
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)
library(MASS)

nperms=100

tis <- args[1]
phen <- args[2]
permset <- as.integer(args[3])

#for testing:
#tis <- "DGN-WB"
#phen <- "ordCIPN8"
#permset <- 0

my.dir <- "/group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/"
#my.dir <- "/Volumes/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/"
exp.dir <- my.dir %&% "genotypes/UMich_imputation_results/PrediXcan/pred_exp/"
pheno.dir <- my.dir %&% "genotypes/"
res.dir <- my.dir %&% "genotypes/UMich_imputation_results/PrediXcan/PrediXcan_results/"

expfile <- exp.dir %&% "N88_n953_" %&% tis %&% ".txt"
phenofile <- pheno.dir %&% "PtStudy.phenofile"
covfile <- pheno.dir %&% "PtStudy.covfile"
gencodefile <- exp.dir %&% "gencode.v18.genes.patched_contigs.summary.protein"

########functions
##adapted from http://www.ats.ucla.edu/stat/r/dae/ologit.htm
ordreg <- function(gt,pt=phenotype,covariates=covmat){
  m<-polr(pt ~ gt + covariates, Hess=TRUE)
  p <- pnorm(abs(summary(m)$coefficients[, "t value"]), lower.tail = FALSE) * 2
  or <- exp(summary(m)$coefficients[,"Value"])
  ci <- exp(confint.default(m))
  colnames(ci) <- c("lCI","uCI")
  n <- m$n
  names(n) <- "n"
  res <- cbind(summary(m)$coefficients[,],"OR"=or,"pval"=p)[1,]
  res <- c(res, ci[1,], n)
  return(res)
}
#just return pval for gene
shortordreg <- function(gt,pt=phenotype,cov=covmat){
  m<-polr(pt ~ gt + cov, Hess=TRUE)
  p <- pnorm(abs(summary(m)$coefficients[, "t value"]), lower.tail = FALSE) * 2
  return(p[1])
}

exp <- read.table(expfile,header=TRUE)
###to prevent "fitted probabilities numerically 0 or 1 occur", the MAF must be > 0.05
#goodsnps <- as.character(readRDS(geno.dir %&% "N88.imputed_MAF0.05_INFOlt1.05_SNPlist.rds"))
#geno <- geno[rownames(geno) %in% goodsnps,]
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

phenotype <- dplyr::select(pheno,starts_with(phen))
phenotype <- factor(phenotype[,1])

gencode <- read.table(gencodefile)
colnames(gencode) <- c("chr","strand","start","end","ensid","gene","type","status")

perm.mat <- matrix(NA,nrow=dim(texp)[1],ncol=nperms)
##run ordinal regression:

###1. covariates: agediagnosis (only agediagnosis associates with ordCIPN8)
covmat <- dplyr::select(pheno,agediagnosis)
covmat <- as.matrix(covmat)
ordgwas <- apply(texp,1,ordreg)
if(tis == "DGN-WB"){
  tordgwas <- data.frame(t(ordgwas)) %>% dplyr::mutate(gene=colnames(ordgwas)) %>% dplyr::mutate(P=pval)
  output <- left_join(tordgwas,gencode,by="gene") %>% dplyr::select(-type,-status,-pval) 
  output <- output[duplicated(output$gene)==FALSE,]
  }else{
  tordgwas <- data.frame(t(ordgwas)) %>% dplyr::mutate(ensid=colnames(ordgwas)) %>% dplyr::mutate(P=pval)
  output <- left_join(tordgwas,gencode,by="ensid") %>% dplyr::select(-type,-status,-pval) 
}

seed = permset*42
data = cbind(phenotype,covmat)
for(i in 1:nperms){
  d <- data[base::sample(1:length(phenotype)),]
  phenotype <- factor(d[,1])
  covmat <- as.matrix(d[,2:dim(d)[2]])
  permpx <- apply(texp,1,shortordreg)
  perm.mat[,i] <- signif(permpx,4)
}

#calc emp P 
a<-output$P >= perm.mat #>= counts zeros
b<-output$P > perm.mat #> doesn't count zeros
empPa <- rowSums(a)/dim(a)[2]
empPb <- rowSums(b)/dim(b)[2]
output <- mutate(output, empPa=empPa, counta=rowSums(a), empPb=empPb, countb=rowSums(b))

write.table(output, file=res.dir %&% "N88_n953_" %&% phen %&% "_agediagnosis_" %&% tis %&% ".ordreg.PrediXcan.empP_permset" %&% 
              permset,quote=F,row.names=F)
