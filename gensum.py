import sys
import os

def cat(suffix, oname, ftype):
  with open(oname+'.'+ftype, 'w') as outfile:
    firstfile = True
    for fname in os.listdir('results/'+ftype):
      if fname.startswith(suffix):
        with open('results/'+ftype+'/'+fname) as infile:
          if (firstfile == False):
            infile.readline()
          outfile.write(infile.read())
        firstfile = False


def main():
  args    = sys.argv[1:]
  suffix  = args[0]
  oname   = args[1]
  cat(suffix, oname, 'char')
  cat(suffix, oname, 'pvalues')


if __name__ == '__main__':
  main()