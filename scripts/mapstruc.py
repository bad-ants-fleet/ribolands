#!/usr/bin/env python

""" 
"""

import os
import re
import sys
import argparse

import ribolands as ribo


def add_sys_barriers_args(parser):
  parser.add_argument("--barriers", default='barriers', action = 'store',
      help="Specify path to your *barriers* executable") 

  parser.add_argument("-e", "--s_ener", type=float, default=None,
      help="Set the energy range for suboptimal structure computation." + 
      " This will overwrite the options --s_maxe and --s_maxn.")
  parser.add_argument("--s_maxe", type=float, default=20,
      help="Set a the maximum subopt range in combination with --s_maxn.")
  parser.add_argument("--s_maxn", type=int, default=7100000,
      help="Specify the number of suboptimal structures. The corresponding" +
      " energy range is computed form the full length molecule.")
  parser.add_argument("--b_minh", type=float, default=0.001,
      help="Set the minimum barrier height (i.e. barriers --minh)")
  parser.add_argument("--b_maxn", type=int, default=100,
      help="Set the maximum number of local minima (i.e. barriers --max)")

def scalar_product(pX, pY, cutoff, verbose=False):
  """calculate the similarity score of two simulations"""
  nX, nY = [], []
  for lmin in range(1, len(pX)):
    if pX[lmin] > cutoff or pY[lmin] > cutoff :
      nX.append(pX[lmin])
      nY.append(pY[lmin])

  nX = rescale(nX, 1.0)
  nY = rescale(nY, 1.0)

  if args.verbose :
    print "pX:", ' '.join(["{:.5f}".format(float(x)) for x in pX]) 
    print "pX:", ' '.join(["{:.5f}".format(float(x)) for x in nX]) 
    print "pY:", ' '.join(["{:.5f}".format(float(y)) for y in pY]) 
    print "pY:", ' '.join(["{:.5f}".format(float(y)) for y in nY]) 

  return sum([sqrt(x*y) for x,y in zip(pX,pY)])

def barriersCG(args, name, seq) :
  mapper = dict() 
  with open(args.mfile) as m :
    for e, line in enumerate(m) :
      if re.match('[\.\(\)]+', line.strip()):
        data = line.split()[:2]
        ss = data[0]
        oc = float(data[1]) if data[1] else 0.
        mapper[e]=[ss,oc]
      else :
        raise ValueError("mapfile does not have the required format.")

  sfile = ribo.sys_suboptimals(name, seq, 
      ener=args.s_ener, 
      noLP=args.noLP,
      temp=args.temperature,
      verb=args.verbose, 
      force=args.force)

  [sfile, bfile, efile, rfile, psfile] = ribo.sys_barriers(name, seq, sfile, 
      barriers=args.barriers,
      minh=args.b_minh,
      maxn=args.b_maxn,
      temp=args.temperature,
      noLP=args.noLP,
      moves='single-base-pair',
      gzip=True,
      rates=True,
      bsize=False,
      saddle=False,
      mfile=args.mfile,
      force=True,
      verb=args.verbose)

#  print coarsify(bfile, 3) 

  return get_pop_vector(mapper, efile, bfile)

#def coarsify(bfile, minh):

 
def get_pop_vector(id_oc, mfile, bfile):
  """Return a vector of local minima and their occupancy. """

  lm = [] # lm[idx] = [structure, occupancy]
  with open(bfile) as tree:
    for n, line in enumerate(tree):
      if n == 0 : # Sequence
        lm.append([None, None])
      else : # old idx = enumerate n
        lm.append([line.split()[1], 0.])

  with open(mfile, 'r') as mapf:
    n = 0
    for line in mapf:
      if re.match('[\.\(\)]+', line.strip()):
        gstr, sptidx, energy, fmin, fminT, gmin, gminT = line.strip().split()
        oc = id_oc[n][1]
        lm[int(gminT)][1] += oc
        del id_oc[n]
        n += 1
      elif re.match('not in hash', line.strip()):
        #raise ValueError("structure not in hash")
        print "WARNING: structure not in hash"
        n += 1
      elif re.match('not yet assigned', line.strip()):
        print "WARNING: structure not yet assigned"
        n += 1

  if sum(id_oc.values()) :
    raise Exception("remove this line <-")
    lm[0][1] = sum(id_oc.values())
  return lm


def main():
  """Compute a population vector using barrier-tree coarse-graining.

  A wrapper that reads a sequence and a list of secondary structures and their
  occupancy to map them into the coarse-graining of a barrier-tree. The
  occupancy of all specified stats should sum up to 1. `mapstruc.py` will
  return a vector of leaves of the corresponding barrier-tree, and the
  occupancy of each basin.

  Input: Sequence or fasta-file.
  Required argument: --mfile (Format = Structure Occupancy ...)
    ((((((...(((....)))))))))(((((((..........))))))).  0.870 # -7.80    0   3.00
    ((((((.............))))))(((((((..........))))))).  0.100 # -6.50    2   1.10
    .(((((.(.(((((....))))).))))))....................  0.010 # -5.90    0   1.10
    ((((((...(((....))))))))).........................  0.005 # -5.60    0   0.80
    ((((((..((.....))..)))))).........................  0.010 # -5.20    0   0.40
    ....((((.(((((....))))).))))((((..........))))....  0.005 # -5.20    0   0.40

  Optional arguments: Arguments for barriers, ViennaRNA, etc.

  """
  # Example for input:
  # randseq -l 40 | DrTransformer.py --k0 1 --t8 4000 -tX 4000 | grep LAST | awk '{print $3, $5}'

  ### parse input & adjust arguments ###
  parser = argparse.ArgumentParser()
  parser.add_argument("--mfile", required=True, action = 'store')

  parser.add_argument("--tmpdir", default='BarMap_Files', action = 'store',
      help="Specify path for storing temporary output files")
  parser.add_argument("-v","--verbose", action="store_true",
      help="Track process by writing verbose output during calculations") 
  parser.add_argument("-f","--force", action="store_true",
      help="Force to overwrite existing BarMap files") 
  parser.add_argument("-n","--name", default='',
      help="Name your output files, this option overwrites the fasta-header")

  parser.add_argument("--noLP", action="store_true",
      help="Compute BarMap pipeline with the ViennaRNA --noLP option")
  parser.add_argument("-T","--temperature", type=float, default=37.0,
      help="Set the temperature for ViennaRNA computations")

  #add_sys_suboptimals_args(parser)
  add_sys_barriers_args(parser)
  args = parser.parse_args()

  name, seq = ribo.parse_vienna_stdin(sys.stdin)

  lm = barriersCG(args, name, seq)

  print "Population Vector with Barriers-mapping"
  print sum(map(lambda x:x[1], lm[1:])), map(lambda x:x[1], lm)
  if args.verbose:
    for l in lm: print l
  return

if __name__ == '__main__':
  main()


