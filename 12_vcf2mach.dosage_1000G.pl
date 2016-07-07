#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw[min max];
use Compress::Zlib;

####This perl script takes the UMich imputed vcf files (MAF>0.01) as input, 
#### removes SNPs with R2<0.8, finds the rsID for each SNP, and makes several output files for each autosome for future parallel computing:
#### .mlinfo.gz and .mldose.gz MACH files for GCTA
#### .dose.gz for plink 
#### .SNPxID.rds for quick reading into R
#### .bim plink bim file with SNP pos info in .dose.gz #dose allele is A1


my $dir = "/group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/";
my $refdir = "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/DGN-WB_genotypes/";

#parse by chr
for(my $i = 1; $i <= 22; $i++){
    print "$i\n";
    my $snpxidhandle = "SNPxID" . $i;
    my $mlinfohandle = "MLINFO" . $i;
    my $bimhandle = "BIM" . $i;
    open($snpxidhandle, ">${dir}mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID");
    open($bimhandle, ">${dir}mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim");
    open($mlinfohandle, ">${dir}mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.mlinfo");
    open(ID, ">${dir}mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.ID.list");
    
    open(REF, "${refdir}ALL.chr${i}.SHAPEIT2_integrated_phase1_v3.20101123.snps_indels_svs.all.noSingleton.RECORDS");

    my %rsid;
    while(<REF>){
	chomp;
	my ($chrom, $pos, $rs) = split(/\t/);
	$rsid{$pos} = $rs;
    }

    open(INTRO, ">intro");

    my $gz = gzopen("${dir}chr${i}.dose_maf0.01_rm.indel.vcf.gz","rb") or die "Error reading vcf.gz file\n";
    while($gz->gzreadline(my $line) > 0){
	chomp($line);
	if ($line =~ /^##/) {next;} # skip vcf header
	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genos) = split(/\t/,$line);
	my ($expfreq, $impr2, $imper2) = split(/;/,$info);
	my ($b, $rsq) = split(/=/,$impr2);
	my $quality = $rsq; ##only one quality score per SNP in minimac, output twice in mach files 
	my $rs = $rsid{$pos}; ##get rsid from .RECORDS hash

	if($chr eq "#CHROM"){
	    my $mlinfohandle = "MLINFO" . $i;
	    print $mlinfohandle "SNP\tAl1\tAl2\tFreq1\tMAF\tQuality\tRsq\n";
	
	    foreach my $id (@genos){
		my (@ids) = split(/_/,$id);
		my $shortid = $ids[0] . "_" . $ids[1] . "_" . $ids[2]; #specific to N88 samples
		print ID "$shortid\n";
		print INTRO "$shortid->$shortid MLDOSE\n";
	    }
	}elsif($rs =~ m/rs\d+/ && $rsq >= 0.8 && $pos =~ m/\d+/){ ###only pull SNPs with rsIDs & R2>0.8 & don't print header rows
	    my $snpxidhandle = "SNPxID" . $chr;		
	    my $bimhandle = "BIM" . $chr;
	    my $mlinfohandle = "MLINFO" . $chr;
	    my $sum = 0; #for $freqalt calc
	    foreach my $geno (@genos){
		my ($gt, $dos) = split(/:/,$geno);
		print $snpxidhandle "$dos\t";
		$sum = $sum + $dos; #add up the dosages to calc ALT allele freq
	    }
	    my $n = @genos; #n samples
	    my $freqalt = $sum/($n*2); #calc ALT allele freq (I found that ALT is not always the minor allele)
	    my $freqref = 1 - $freqalt;
	    my $maf = min($freqref,$freqalt);
	    $maf = sprintf("%.5f", $maf); #round to 5 places
	    $freqref = sprintf("%.5f", $freqref);
	    print $bimhandle "$chr\t$rs\t0\t$pos\t$alt\t$ref\n"; #plink recognizes dosages of the A1 allele (ALT allele in vcf)
            print $mlinfohandle "$rs\t$ref\t$alt\t$freqref\t$maf\t$quality\t$rsq\n";
	    print $snpxidhandle "\n";
	}
    }
}


for(my $i = 1; $i <= 22; $i++){
    print "$i\n";
    open(R, ">runR.R") or die "cant make runR.R\n";
    print R "bim<-read.table(\"${dir}mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim\")\n";
    print R "saveRDS(bim,\"${dir}mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim.rds\")\n";
    print R "dat<-scan(\"${dir}mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID\")\n";
    print R "gtidlist<-scan(\"${dir}mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.ID.list\",\"character\")\n";
    print R "dat<-matrix(dat, ncol=length(gtidlist), byrow=T)\n";
    print R "colnames(dat) <- gtidlist\n";
    print R "rownames(dat) <- bim[,2]\n";
    print R "saveRDS(dat,\"${dir}mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.rds\")\n";
    print R "dose <- cbind(bim[,2],bim[,5:6],dat)\n";
    print R "write.table(dose,\"${dir}mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose\",quote=F,row.names=F,col.names=F)\n";

    print R "dat<-t(dat)\n";
    print R "write.table(dat,\"t.dos.chr${i}\",col.names=F,row.names=F,quote=F)\n";
    close(R);
    system("R --vanilla < runR.R");
    system("paste -d\' \' intro t.dos.chr${i} > ${dir}mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.mldose");
    system("cut -f 2,5-6 ${dir}mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.bim >o");
    system("paste o ${dir}mach_dosage_files/N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID > N88.imputed_maf0.01_R20.8_1000G.chr${i}.SNPxID.dose");
}


system("gzip ${dir}mach_dosage_files/*.mldose");
system("gzip ${dir}mach_dosage_files/*.mlinfo");
system("gzip ${dir}mach_dosage_files/*.SNPxID.dose");

#do this after verify everything worked correctly
#system("rm intro t.dos.chr* runR.R ${dir}mach_dosage_files/*1000G.chr*.SNPxID");

