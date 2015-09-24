import sys
import os.path

# given r, create haplotype for rth region
def genhap(r):
  
  regions = open('/n/home01/godwin/Research/ipwskat/data/cosi/10Kregions.txt','r')
  i = 1
  for rline in regions:
    
    if (i == r):
      
      filename = '/n/home01/godwin/Research/ipwskat/data/cosi/10Khaplotypes/haplotypes'+str(r)+'.txt'
      if (not os.path.isfile(filename)):
        # get min and max index for variants within the chosen subregion
        [mini, maxi] = rline.rstrip().split()
        [mini, maxi] = map(int, [mini, maxi])
        
        # create haplotype matrix
        haplotypes = open(filename,'w')
        cosi = open('/n/home01/godwin/Research/ipwskat/data/cosi/out.hap-1','r')
        for cline in cosi:
          hap = cline.rstrip().split()[mini+1:maxi+2] # restrict each line to only variants within the chosen subregion
          for item in hap:
            haplotypes.write(str(2-int(item))+'\t')
          
          haplotypes.write('\n')
        
        haplotypes.close()
      
    i = i + 1


def main():
  args = sys.argv[1:]
  
  taskID1   = int(args[0])
  taskID2   = int(args[1])
  njobs     = int(args[2])
  interval  = 1000/njobs
  rmin      = 10000/njobs*(taskID1-1)+interval*(taskID2-1)+1
  rmax      = 10000/njobs*(taskID1-1)+interval*(taskID2-1)+interval

  file = open('../output/done.txt','r')
  done = file.readline().split()
  file.close()
  
  for r in range(rmin,rmax+1):
      if (str(r) not in done):
          genhap(r)


if __name__ == '__main__':
  main()
