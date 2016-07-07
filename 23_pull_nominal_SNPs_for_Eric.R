##pull nominal SNPs (P<0.05) for Eric

library(dplyr)
library(tidyr)
library(ggplot2)
"%&%" = function(a,b) paste(a,b,sep="")
my.dir = "/Volumes/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/GWAS_results/"

gwasfile <- my.dir %&% 'N88.imputed_rnGM412_age.cisp.10PCs_chr1-22.assoc.dosage.sorted'
gwas <- read.table(gwasfile,header=T)
gwas <- gwas[complete.cases(gwas),] #rm NAs
gwas <- dplyr::filter(gwas,INFO<1.05,P<0.05)
write.table(gwas,file=gwasfile %&% ".P_0.05",row.names = F,quote=F)

gwasfile <- my.dir %&% 'N88.imputed_dosegroup_10PCs_chr1-22.assoc.dosage.sorted'
gwas <- read.table(gwasfile,header=T)
gwas <- gwas[complete.cases(gwas),] #rm NAs
gwas <- dplyr::filter(gwas,INFO<1.05,P<0.05)
write.table(gwas,file=gwasfile %&% ".P_0.05",row.names = F,quote=F)

gwasfile <- my.dir %&% 'N88.imputed_ord3CIPN8_agediagnosis_chr1-22.ordreg.assoc.dosage.sorted'
gwas <- read.table(gwasfile,header=T)
gwas <- gwas[complete.cases(gwas),] #rm NAs
gwas <- dplyr::filter(gwas,P<0.05)
write.table(gwas,file=gwasfile %&% ".P_0.05",row.names = F,quote=F)

