
import re
import sys

def parse_vienna_stdin():
  """ Read STDIN in fasta format
  Read a Sequence and its Name in Fasta Format. Only one Input-Sequence is
  allowed at a time. The Characters must be A, C, G, U, &

  :return: (name, sequence)
  """
  name = 'NoName'
  seq  = ''
  for line in sys.stdin:
    if re.match('>', line):
      if name != 'NoName' :
        print >> sys.stderr, 'Only single-sequence fasta format supported!'
        raise SystemExit
      else : 
        name = line.strip().split()[0][1:]
    else:
      seq += line.strip()
  m = re.search('[^AUCG&]', seq) 
  if m :
    print >> sys.stderr, \
      "Does not look like RNA:", m.string[m.span()[0]], "in", seq
    raise SystemExit
  return (name, seq)

# make sure that you use args in order to name every column correctly
# maybe even return a pandas.DataFrame (?)
def parse_barfile(bfile, seq=''):
  """ return the content of a barriers output-file """
  output = []
  with open(bfile) as bar :
    for e, line in enumerate(bar) :
      if e == 0 : 
        if seq and seq != line.strip() :
          print >> sys.stderr, 'Wrong sequence', seq, ' vs. ', line
          raise SystemExit
      else :
        output.append(line.strip().split())
        #[idx, lmin, en, father, bar] = line.strip().split()[0:5]
        #output.append([idx, lmin, en, father, bar])
  return output

def parse_ratefile(rfile):
  """ return the content of a barriers rates-file """
  RM = []
  with open(rfile) as rates :
    for line in rates :
      RM.append((map(float, line.strip().split())))
  return RM

