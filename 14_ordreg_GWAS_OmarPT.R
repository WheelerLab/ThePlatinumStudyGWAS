########################################
#### Ordinal regression GWAS        ####
#### by Heather E. Wheeler 20151207 ####
########################################
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)
library(MASS)

chr <- args[1]
phen <- args[2]

#for testing:
#chr <- "22"
#phen <- "ordCIPN8"

my.dir <- "/group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/"
#my.dir <- "/Volumes/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/"
geno.dir <- my.dir %&% "genotypes/UMich_imputation_results/mach_dosage_files/"
pheno.dir <- my.dir %&% "genotypes/"
res.dir <- my.dir %&% "GWAS_results/"

genofile <- geno.dir %&% "N88.imputed_maf0.01_R20.8_1000G.chr" %&% chr %&% ".SNPxID.rds"
phenofile <- pheno.dir %&% "CIPN_phenotypes_covariates_fromOmar.txt"
bimfile <- geno.dir %&% "N88.imputed_maf0.01_R20.8_1000G.chr" %&% chr %&% ".bim.rds"

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

geno <- readRDS(genofile)
###to prevent "fitted probabilities numerically 0 or 1 occur", the MAF must be > 0.05
goodsnps <- as.character(readRDS(geno.dir %&% "N88.imputed_MAF0.05_INFOlt1.05_SNPlist.rds"))
geno <- geno[rownames(geno) %in% goodsnps,]
#for testing:
#geno <- geno[1:10,]

bim <- readRDS(bimfile)
colnames(bim) <- c("CHR","SNP","CM","BP","A1","A2")
bim <- mutate(bim, SNP=as.character(SNP))
bim <- bim[bim[,2] %in% goodsnps,]

pheno <- read.table(phenofile,header=T)

phenotype <- dplyr::select(pheno,starts_with(phen))
phenotype <- factor(phenotype[,1])

##run ordinal regression:

###1. covariates: age 
#covmat <- dplyr::select(pheno,age)
#covmat <- as.matrix(covmat)
#ordgwas <- apply(geno,1,ordreg)
#tordgwas <- data.frame(t(ordgwas)) %>% dplyr::mutate(SNP=colnames(ordgwas)) %>% dplyr::mutate(P=pval)

#output <- left_join(bim,tordgwas,by="SNP") %>% dplyr::select(-CM,-pval)

#write.table(output, file=res.dir %&% "N88.imputed_" %&% phen %&% "_age_chr" %&% chr %&% ".ordreg.assoc.dosage",quote=F,row.names=F)

##2. covariates: age, CumlCispdose, PCs1-10
#covmat <- dplyr::select(pheno,age,CumlCispdose,starts_with("PC"))
#covmat <- as.matrix(covmat)
#ordgwas <- apply(geno,1,ordreg)
#tordgwas <- data.frame(t(ordgwas)) %>% dplyr::mutate(SNP=colnames(ordgwas)) %>% dplyr::mutate(P=pval)

#output <- left_join(bim,tordgwas,by="SNP") %>% dplyr::select(-CM,-pval)

#write.table(output, file=res.dir %&% "N88.imputed_" %&% phen %&% "_age.cisp.10PCs_chr" %&% chr %&% ".ordreg.assoc.dosage",quote=F,row.names=F)

##3. covariates: age ever_smoke hbpmed
covmat <- dplyr::select(pheno,age,ever_smoke,hbpmed)
covmat <- as.matrix(covmat)
ordgwas <- apply(geno,1,ordreg)
tordgwas <- data.frame(t(ordgwas)) %>% dplyr::mutate(SNP=colnames(ordgwas)) %>% dplyr::mutate(P=pval)

output <- left_join(bim,tordgwas,by="SNP") %>% dplyr::select(-CM,-pval)

write.table(output, file=res.dir %&% "N88.imputed_" %&% phen %&% "_age.eversmoke.hbpmed_chr" %&% chr %&% ".ordreg.assoc.dosage",quote=F,row.names=F)
