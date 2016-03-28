#!/usr/bin/env python

#  Coded by: Stefan Badelt <stef@tbi.univie.ac.at>
#  University of Vienna, Department of Theoretical Chemistry
#
#  -*- Style -*- 
#  Use double quotes or '#' for comments, such that single quotes are available
#  for uncommenting large parts during testing
#
#  *) do not exceed 80 characters per line
#  *) indents: 2x whitespace, no tab characters!
#
#  -*- VIM config -*- 

#  set textwidth=80
#  set ts=2 et sw=2 sts=2
#
#  -*- Content -*-
#  *) parsers for stdin, barfiles and rate-matrix
#
#  -*- TODO -*-
#  *) write documentation
#  *) add general utilities: make_pair_table

import re
import sys

def make_pair_table(ss, base=0):
  """ 
    Return a secondary struture in form of pair table:
     base=0: ((..)). => [5,4,-1,-1,1,0,-1]
      i.e. start counting from 0, unpaired = -1
     base=1: ((..)). => [7,6,5,0,0,2,1,0]
      i.e. start counting from 1, unpaired = 0, pt[0]=len(ss)

    TODO: raise error for unbalanced brackets
  """
  stack=[];

  if base is 0:
    pt=[-1] * len(ss);
  else :
    base = 1
    pt = [0] * (len(ss) + base);
    pt[0] = len(ss);

  for i, char in enumerate(ss, base):
    if (char == '('):
      stack.append(i);
    elif (char == ')'):
      j=stack.pop();
      pt[i]=j
      pt[j]=i
  return pt

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
        raise ValueError
      else : 
        name = line.strip().split()[0][1:]
    else:
      seq += line.strip()
  m = re.search('[^AUCG&]', seq) 
  if m :
    print >> sys.stderr, \
      "Does not look like RNA:", m.string[m.span()[0]], "in", seq
    raise ValueError
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
          raise ValueError
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

