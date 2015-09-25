import sys
import os.path

# Let S={{1,...,10},{11,...20},...}. Given r, return the index of the subset in S containing r. 
def mapS(r):
  return (r-1)/10+1

# given r, create haplotype for rth region
def genhap(infile1, infile2, folder, min10r, max10r):
  
  regions = open(infile2,'r')
  r = 1
  for rline in regions:
    outfile = folder + '/haplotypes' + str(r) + '.txt'
    if (mapS(r) >= min10r and not os.path.isfile(outfile)):
      
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
    if (mapS(r) > max10r):
      break
  
  regions.close()


def main():
  args = sys.argv[1:]
  filename_hap = args[0] # file containing haplotypes across entire 1Mb region
  filename_reg = args[1] # file containing chromosomal regions (beginning and end bp) to restrict the haplotypes to
  foldname_out = args[2] # folder where restricted haplotypes will be saved to
  if len(args) == 5:
    min10r = int(args[3])
    max10r = int(args[4])
  else:
    min10r = 1
    max10r = float("inf")
  genhap(filename_hap, filename_reg, foldname_out, min10r, max10r)


if __name__ == '__main__':
  main()