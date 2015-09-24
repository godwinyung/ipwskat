# ipwskat

_genhap.py_
Usage: python genhap.py infile1 infile2 outfolder [min] [max]
Description: Given infile1 (rows of haplotypes) and infile2 (rows of region range), for the rth region (min <= r <= max), restrict each haplotype to the region and save the rows of restricted haplotypes to outfolder. If min and max are not specified, then program generates restricted haplotypes for every region.
