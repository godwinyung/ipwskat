import sys
import os.path

# given r, create haplotype for rth region
def genhap(infile1, infile2, folder, minr, maxr):
  
  regions = open(infile2,'r')
  r = 1
  for rline in regions:
    outfile = folder + '/haplotypes' + str(r) + '.txt'
    if (r >= minr and not os.path.isfile(outfile)):
      
      # get min and max index for variants within the chosen subregion
      [mini, maxi] = rline.rstrip().split()
      [mini, maxi] = map(int, [mini, maxi])
      
      # create haplotype matrix
      out = open(outfile,'w')
      haplotypes = open(infile1,'r')
      for haplotype in haplotypes:
          rhap = haplotype.rstrip().split()[mini+1:maxi+2] # restrict each line (haplotype) to only variants within the chosen subregion
          for item in rhap:
            out.write(str(2-int(item))+'\t')
          out.write('\n')
      out.close()
      haplotypes.close()
      
    r = r + 1
    if (r > maxr):
      break
  
  regions.close()


def main():
  args = sys.argv[1:]
  filename_hap = args[0] # file containing haplotypes across entire 1Mb region
  filename_reg = args[1] # file containing chromosomal regions (beginning and end bp) to restrict the haplotypes to
  foldname_out = args[2] # folder where restricted haplotypes will be saved to
  if len(args) == 5:
    minr = int(args[3])
    maxr = int(args[4])
  else:
    minr = 1
    maxr = float("inf")
  genhap(filename_hap, filename_reg, foldname_out, minr, maxr)


if __name__ == '__main__':
  main()
