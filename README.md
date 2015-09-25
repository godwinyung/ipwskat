# ipwskat

**genhap.py**  
_Usage_: python genhap.py infile1 infile2 outfolder [min] [max]  
_Description_: Let S={{1,...,10},{11,...,20},...} be the set of 10 consecutive natural numbers. Given infile1 (rows of haplotypes) and infile2 (rows of region range), for regions {S[min+i][j]} (i=0,...,max-min and j=1,...,10), restrict each haplotype to the region and save the rows of restricted haplotypes to outfolder. If min and max are not specified, then restricted haplotypes are generated and saved for every region, i.e. min=1, max=inf.

**genhap.sbatch**

**ipwskat.R**

**sim-cts.R**

**sim-cts.sbatch**

**sim-setting.R**
