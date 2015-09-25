# ipwskat

**genhap.py**  
_Usage_: python genhap.py infile1 infile2 outfolder [min] [max]  
_Description_: Let S={{1,...,10},{11,...,20},...} be the set of 10 consecutive natural numbers. Given infile1 (rows of haplotypes) and infile2 (rows of region range), for regions {S[min+i][j]} (i=0,...,max-min and j=1,...,10), restrict each haplotype to the region and save the rows of restricted haplotypes to outfolder. If min and max are not specified, then restricted haplotypes are generated and saved for every region, i.e. min=1, max=inf.

**genhap.sbatch**
_Usage_: sbatch --array=1-1000 genhap.sbatch
_Description_: Generates 10,000 haplotypes for each of 10,000 regions.

**ipwskat.R**
_Usage_: IPWSKAT_Null_Model(), IPWSKAT()
_Description_: See documentation for SKAT.

**sim-cts.R**
_Usage_: sim.cts()
_Description_: Given a set of parameters, SNP information, haplotypes, and output file name, simulates data under the parameter settings and saves results as an R object with:
1. table of p-values
2. table of population statistics (e.g., prevalence, number of observed variants)
Note, results are saved intermittently so that jobs can be stopped and restarted at any time.

**sim-cts.sbatch**

**sim-setting.R**
