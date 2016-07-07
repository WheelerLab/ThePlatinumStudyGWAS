#!/usr/bin/env python

file =  open("/group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/N88_Recluster_TOP_20150911_FinalReport.postPCA.euro.fam")
outfile = open("/group/dolan-lab/hwheeler/ThePlatinumStudy/GWAS/genotypes/plate_check_phenotypes","w")

outfile.write("FID IID plate1 plate2 plate3 plate4 plate5 plate6 plate7 plate8 plate9 plate10 plate11 plate12\n")

for line in file:
    data = line.split(' ')
    outfile.write(data[0] + ' ' + data[1])
    plate = data[0].split('_')[1]
    for i in range(1,13):
        platenum = "Plate" + str(i).zfill(2)
        if(plate == platenum):
            outfile.write(' ' + str(2))
        else:
            outfile.write(' ' + str(1))
    outfile.write('\n')


