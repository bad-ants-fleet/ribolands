#!/usr/bin/env python

""" 
  Written by Stefan Badelt (stef@tbi.univie.ac.at)

  Department of Theoretical Chemistry, University of Vienna
  http://www.tbi.univie.ac.at

  vim-config = set: ts=2 et sw=2 sts=2
"""

import re
import sys
import argparse
import string
import math
import RNA

import ribolands as ril

def aptamer_energy(seq, ss, verb=False,
    # Default Theophylline
    apt = 'GAUACCAG'+'&'+'CCCUUGGCAGC',
    poc = '(...(((('+'&'+')...)))...)',
    bfe = -8.86): # at 25*C; -9.22 at 37*C
  """ 
    Check if a sequence/structure pair contains the 
    ligand binding pocket (apt/poc). If so, return the 
    binding free energy (bfe), otherwise return 0.
    Multiple pockets will return bfe multiple times!

    TODO: allow hairpin pockets (e.g. tetracycline)
  """

  [aptL, aptR] = apt.split('&')
  [pocL, pocR] = poc.split('&')

  patL = re.compile(aptL)
  patR = re.compile(aptR)

  sites=0
  for mL in patL.finditer(seq) :
    if pocL != ss[mL.start():mL.end()] :
      continue
    for mR in patR.finditer(seq) :
      if mR.start() < mL.end() :
        continue
      if pocR != ss[mR.start():mR.end()] :
        continue
      # Now make sure that these are really base-pairs
      ptable = ril.make_pair_table(ss, base=0)
      if mL.start() == ptable[mR.end()-1] and \
          mL.end()-1 == ptable[mR.start()] and \
          mR.start() == ptable[mL.end()-1] and \
          mR.end()-1 == ptable[mL.start()] :
        #if verb :
        #  print >> sys.stderr, "{:s} {:6.2f}".format(ss, bfe)
        sites += 1
  return bfe * sites

def check_symmetry(sym, seq) :
  """ Check if a cofolded sequence is a homo-dimer """
  if seq.find('&') == -1 :
    return False
  [s1,s2] = seq.split('&')
  if sym :
    return (s1 == s2)
  elif (s1 == s2) :
    print >> sys.stderr, "spatch.py: Ignoring symmetry correction for homo-dimer (see Option -s)!"
  return False

def is_symmetric(ss):
  """ See if a homo-dimer secondary structure has rotational symmetry
    Test e.g. with "CCCGGCCGGG&CCCGGCCGGG" suboptimals

    Of the following secondary structures:
    1 .(.(..(.(.&.).)..).).
    2 .(.(..(.).&.(.)..).).
    3 .(.(..(.).&.).)..(.).
    4 .(.)..(.).&.(.)..(.).

    5 .(.(..(.)..(.(.&.).)..(.)..).).
    6 .(.(..(.).(.(.&.).)..(.).).).
    7 .(.(..(.).(.(.&.).).(.)..).).

    only 1 and 5 are considered as rotational symmetric, 
    only 6 and 7 are filtered by the fast exclusion of rotational symmetry

    :return: True/False
  """
  [s1,s2] = ss.split('&')

  """ fast way to exclude rotational symmetry """
  t1 = s1.translate(string.maketrans("()", "xx"))
  t2 = s2.translate(string.maketrans("()", "xx"))
  if t1 != t2 or t1 != t1[::-1] : 
    return False

  """ slow way to ensure rotational symmetry """
  stack = []
  pt=[0] * len(s1);
  for j, char in enumerate(s1) :
    if char == '(' : 
      stack.append(j)
    elif char == ')' : 
      i = stack.pop()
      pt[i]=j
      pt[j]=i

  if not stack : 
    """ Got a fake-dimer, no '&' crossing base-pairs """
    return False
  else :
    for i in stack : pt[i] = None

  palin = ''
  for i in pt :
    if i == None : 
      palin += '|'
    elif i == 0 : 
      palin += '.'
    else :
      palin += 'x'
  return (palin == palin[::-1])

def cofold_noLP_energy(ss, seq, noLC=False):
  """ Correct cofold energies returned from RNAsubopt """
  if noLC and ((re.search(re.escape('.(&'), ss) or \
      re.search(re.escape('&).'), ss) or \
      re.search(re.escape('&(.'), ss) or \
      re.search(re.escape('.)&'), ss))) :
    return None
  if seq :
    bup = RNA.cvar.cut_point
    [r1,r2]=seq.split('&')
    [s1,s2]=ss.split('&')
    RNA.cvar.cut_point = len(r1)+1
    en = RNA.energy_of_structure(r1+r2, s1+s2, 0)
    RNA.cvar.cut_point = bup
  return en

def is_true_dimer(ss):
  """ return True if there is at least one 
  base pair crossing the '&' character 
  """
  o, c = 0, 0
  for char in ss :
    if char == '(' : o += 1
    elif char == ')' : c += 1
    elif char == '&' : break
  return o != c

def main():
  """ A collection of utils to modify the output of RNAsubopt 
    TODO: tetracycline binding pockets
    TODO: logarithmic multiloops?
    TODO: add temperature
  """
  parser = argparse.ArgumentParser()
  parser.add_argument("--theophylline",
    help="Add an energy term (-8.86 kcal/mol) to strucutes with a \
    theophylline binding pocket", 
    action="store_true")
  parser.add_argument("-d","--dimers",
    help="Chose to print only true dimers, i.e structures with at least one \
    base pair crossing the '&' character", 
    action="store_true")
  parser.add_argument("-s","--symmetry",
    help="Add an entropic symmetry correction penalty to symmetric homo-dimers \
    to compensate for their two-fold rotational symmetry", 
    action="store_true")
  parser.add_argument("--fix_cofold_noLP",
    help="Remove lonely pairs around the '&' character and correct the energy", 
    action="store_true")
  """ Sthg to consider in the future
  parser.add_argument("-T","--temperature",
    help="Set the temperature for symmetry correction", 
    type=float,
    default=37.0)
  """
  parser.add_argument("-v","--verbose",
    help="Verbose output",
    action="store_true")
  args = parser.parse_args()

  RT=0.61632077549999997
  seq=''
  """ Main loop, parse and modify RNAsubopt """
  for e, line in enumerate(sys.stdin):

    if e == 0 : 
      [seq, mfe, enr] = line.strip().split()
      args.symmetry = check_symmetry(args.symmetry, seq)
      print line.strip()
      continue

    """ Read Structure and Energy """
    [ss, en] = line.strip().split()

    if args.theophylline :
      en = float(en) + aptamer_energy(seq, ss, verb=args.verbose, bfe=-9.32)

    if args.fix_cofold_noLP :
      en = cofold_noLP_energy(ss, seq, True) 
      if en == None :
        continue

    if args.dimers and not is_true_dimer(ss) :
      """ Structure is not a real dimer """
      continue

    #if args.symmetry and is_symmetric(ss) : # how it should be
    if args.symmetry and is_true_dimer(ss) : # consistent with RNAcofold
      """ In order to be consistent with partition function calculations:
        the structures need a symmetry correction in case they are true dimers
        is it a problem that they are counted twice if they are not true dimers?
        => in both cases, the correction factor -RT*ln(2) should do it!  
      """
      en = float(en) + RT * math.log(2)

    print "%s %6.2f" % (ss, float(en))

  return

if __name__ == '__main__':
  main()

