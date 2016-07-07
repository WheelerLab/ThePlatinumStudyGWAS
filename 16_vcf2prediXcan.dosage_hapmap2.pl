#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw[min max];
use Compress::Zlib;

####This perl script takes the UMich imputed vcf files (MAF>0.01) as input, 
#### removes SNPs with R2<0.8, removes ambiguous-strand SNPs (A/T and C/G), 
#### removes non-hapmap2 SNPs, finds the rsID for each SNP, and makes output files 
#### for each autosome for PrediXcan:
#### chr${i}.dosage.txt.gz
#### samples.txt
#### dose allele is Allele2, see https://github.com/hakyimlab/PrediXcan/blob/master/Software/HOWTO-beta.md


my $dir = "/group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/UMich_imputation_results/";
my $refdir = "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/DGN-WB_genotypes/";

##predictdb weights currently limited to hapmap2
open(HAP, "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/hapmapSnpsCEU.list"); ##from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hapmapSnpsCEU.txt.gz

my %hapmapsnps;
while(<HAP>){
    chomp;
    my ($snp) = split(/\n/);
    $hapmapsnps{$snp} = 1;
}

#parse by chr
for(my $i = 1; $i <= 22; $i++){
    print "$i\n";
    my $predhandle = "PREDX" . $i;
    open($predhandle, ">${dir}/PrediXcan/dosages/chr${i}.dosage.txt");
    
    open(REF, "${refdir}ALL.chr${i}.SHAPEIT2_integrated_phase1_v3.20101123.snps_indels_svs.all.noSingleton.RECORDS");

    my %rsid;
    while(<REF>){
	chomp;
	my ($chrom, $pos, $rs) = split(/\t/);
	$rsid{$pos} = $rs;
    }

    open(SAMP, ">${dir}/PrediXcan/dosages/samples.txt");

    my $gz = gzopen("${dir}chr${i}.dose_maf0.01_rm.indel.vcf.gz","rb") or die "Error reading vcf.gz file\n";
    while($gz->gzreadline(my $line) > 0){
	chomp($line);
	if ($line =~ /^##/) {next;} # skip vcf header
	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genos) = split(/\t/,$line);
	my ($expfreq, $impr2, $imper2) = split(/;/,$info);
	my ($b, $rsq) = split(/=/,$impr2);
	my $rs = $rsid{$pos}; ##get rsid from .RECORDS hash

	if($chr eq "#CHROM"){
	    foreach my $id (@genos){
		my (@ids) = split(/_/,$id);
		my $shortid = $ids[0] . "_" . $ids[1] . "_" . $ids[2]; #specific to N88 samples
		print SAMP "$shortid\n";
	    }
	}if($ref eq "A" && $alt eq "T"){ ##rm potentially ambiguous strand SNPs
	    next;
	}elsif($ref eq "T" && $alt eq "A"){
	    next;
	}elsif($ref eq "C" && $alt eq "G"){
	    next;
	}elsif($ref eq "G" && $alt eq "C"){
	    next;
	}					
	elsif(defined($hapmapsnps{$rs}) && $rs =~ m/rs\d+/ && $rsq >= 0.8 && $pos =~ m/\d+/){ ###only pull hapmap2 SNPs with rsIDs & R2>0.8 & don't print header rows
	    my $sum = 0; #for $freqalt calc
	    foreach my $geno (@genos){
		my ($gt, $dos) = split(/:/,$geno);
		$sum = $sum + $dos; #add up the dosages to calc ALT allele freq
	    }
	    my $n = @genos; #n samples
	    my $freqalt = $sum/($n*2); #calc ALT allele freq (I found that ALT is not always the minor allele)
	    $freqalt = sprintf("%.5f", $freqalt); #round to 5 places
	    print $predhandle "chr$chr $rs $pos $ref $alt $freqalt";
	    foreach my $geno (@genos){
                my ($gt, $dos) = split(/:/,$geno);
                print $predhandle " $dos";
	    }
	    print $predhandle "\n";
	}
    }
}


system("gzip ${dir}PrediXcan/dosages/*.dosage.txt");




