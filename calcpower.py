import sys
import os.path
import numpy


def read_result(filename):
  file    = open(filename, 'r')
  header  = file.readline().split()
  analyses= header[2:]
  N = len(analyses)
  result  = {analysis:[] for analysis in analyses}
  for fline in file:
    fline = fline.split()[2:]
    for n in range(N):
      if fline[n] != 'NA':
        result[analyses[n]].append(float(fline[n]))
    
  file.close()
  return analyses, result


def calc_power(result, alpha):
  power = {analysis:str(numpy.mean([int(x<alpha) for x in result[analysis]])) for analysis in result.keys()}
  return power


def gen_output(analyses, power, filename):
  ofilename = filename.split('.')[0]+'.txt'
  ofile = open(ofilename ,'w')
  for analysis in analyses:
    ofile.write(analysis+'\t')
    ofile.write('\t'.join(power[analysis]))
    ofile.write('\n')
  ofile.close()


define main():
  args = sys.argv[1:]
  for filename in args:
    analyses, result = read_result(filename)
    power = {analysis:[] for analysis in analyses}
    for alpha in [0.05, 0.01, 0.001, 0.0001, 2.5e-6]:
      temp = calc_power(result,alpha)
      for analysis in analyses:
        power[analysis].append(temp[analysis])
      
    gen_output(analyses, power, filename)
  

if __name__ == '__main__':
  main()