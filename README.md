# ipwskat

* **genhap.py** <br>
_Usage_: python genhap.py infile1 infile2 outfolder [min] [max] <br>
_Description_: Given infile1 (rows of haplotypes) and infile2 (rows of region range), for the rth region (min <= r <= max), restrict each haplotype to the region and save the rows of restricted haplotypes to outfolder. If min and max are not specified, then restricted haplotypes are generated and saved for every region, i.e. min=1, max=inf.

* **genhap.sbatch**
