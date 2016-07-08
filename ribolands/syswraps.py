#!/usr/bin/env python

#  Coded by: Stefan Badelt <stef@tbi.univie.ac.at>
#  University of Vienna, Department of Theoretical Chemistry

#  -*- Content -*-
#  *) systemcalls of RNAsubopt, barriers and treekin 
#  *) tested under linux

#  -*- Style -*- 
#  Use double quotes or '#' for comments, such that single quotes are available
#  for uncommenting large parts during testing
#
#  *) do not exceed 80 characters per line
#  *) indents: 2x whitespace, no tabs!

#  -*- VIM config -*- 
#  set textwidth=80
#  set ts=2 et sw=2 sts=2

#  -*- TODO -*-
#  *) encapsulate syswraps into tmp-directory (bec of barriers)

import os
import re
import sys
import gzip
import math
import subprocess as sub

class ExecError(Exception):
  """Ribolands-error handling

  Attributes:
  msg (str): Human readable string describing the exception.
  """

  def __init__(self, msg):
    """Exception initialization

    Args:
      msg (str): Human readable string describing the exception.
    """
    self.msg = msg

  def __str__(self):
    return "Error: " + self.value

# **************** #
# public functions #
# ................ #

def sys_treekin(name, seq, bfile, rfile,
  treekin = 'treekin',
  verb = False,
  p0 = ['1=0.5', '2=0.5'],
  t0 = 1e-6,
  repl=None,
  useplusI=True,
  ti = 1.02,
  t8 = 1e10,
  force=False):
  """ **Perform a system-call of the program ``treekin``.**

  Printing the results into a file and return the respective filename. This
  wrapper will produce two output files from ``STDIN`` and ``STDERR``,
  respectively. 

  .. note:: This wrapper is written for ``treekin v0.4``. 

  :param name: Name of the sequence used to name the output file
  :param seq: Nucleic acid sequence 
  :param bfile: input filename for ``treekin`` as produced by ``barriers``
  :param rfile: input filename for ``treekin`` to specify a rate-matrix as \
      produced by ``barriers``
  :type name: string
  :type seq: string
  :type bfile: string
  :type rfile: string

  :param p0: A list to specify an initial occupancy vector 
  :type p0: list

  :return: The name of the file containing ``treekin`` results. 
  :rtype: string
  """

  if which(treekin) is None :
    print treekin, "is not executable"
    print """ 
    You need to install *treekin*, which you can download from the 
    ViennaRNA package homepage: http://www.tbi.univie.ac.at/RNA/Treekin/
    
    If you have installed the program, make sure that the path you specified 
    is executable.
    """
    raise RuntimeError('Could not find executable')


  reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')
  # http://www.regular-expressions.info/floatingpoint.html

  tfile = name+'.tkn'
  efile = name+'.err'

  if not force and os.path.exists(tfile):
    if verb : print >> sys.stderr, tfile, "<= Files exist"
    return tfile
    
  treecall = [treekin, '--method', 'I']
  if useplusI : treecall.extend(['--useplusI'])
  treecall.extend(['--tinc', str(ti)])
  treecall.extend(['--t0', str(t0)])
  treecall.extend(('--t8', str(t8)))
  for p in p0:
    treecall.extend(('--p0', p))
  treecall.extend(('-f', rfile))
  
  if verb :
    print >> sys.stderr, ' '.join(treecall), '<', bfile, '2>', efile, '>', tfile

  # Do the simulation (catch treekin errors)
  with open(bfile, 'r') as bar, \
    open(tfile, 'w') as tkn, \
    open(efile, 'w') as err:
      proc = sub.Popen(treecall,stdin=bar,stdout=tkn,stderr=err)
      proc.communicate(None)
      if proc.returncode :
        print >> sys.stderr, \
            "ERROR: process terminated with return code", proc.returncode
        print >> sys.stderr, \
            ' '.join(treecall), '<', bfile, '2>', efile, '>', tfile
        raise RuntimeError

  '''
  # This should be caught by the proc.returncode before ...
  lastlines = sub.check_output(['tail', '-2', tfile]).strip().split("\n")
  if not reg_flt.match(lastlines[0]):
    t_total += t8
    print "{:3d} {:3d} {:f} {:s}".format(l, 1, 1.0, get_structure(l, 1, args))
    if get_structure(l, 2, args):
      print >> sys.stderr, \
          "No valid time-course in {:s}: {:s}".format(tfile, lastlines[0])
      with open(efile, 'r') as err:
        print >> sys.stderr, err.read().strip()
      sys.exit('over and out')
    continue

  time = float(lastlines[0].split()[0])
  iterations = int(lastlines[-1].split()[-1])
  if verb and iterations < iterate :
    print >> sys.stderr, "Equilibrium reached after" + \
    " {:d} of {:d} iterations at time {:f}:".format(iterations, iterate, time)
  '''

  return tfile

def sys_barriers(name, seq, sfile,
    barriers='barriers',
    minh=0.001,
    maxn=50,
    k0=1.0,
    temp=37.0,
    noLP=False,
    moves='single-base-pair',
    gzip=True,
    rates=True,
    bsize=False,
    circ=False,
    saddle=False,
    mfile='',
    force=False,
    verb=False):
  """ **Perform a system-call of the program ``barriers``.**

  The print the results into a file and return the respective filename. This
  wrapper will produce two output files from ``STDIN`` and ``STDERR``,
  respectively. 

  .. note:: This wrapper is written for ``barriers v1.6``, \
      previous implementations do not have the ``--mapstruc`` option.

  :param name: Name of the sequence used to name the output file.
  :param seq: Nucleic acid sequence.
  :param sfile: Input filename for ``barriers`` as produced by ``RNAsubopt``.
  :param force: Overwrite existing files with the same name.
  :param verb: Print the current system-call to ``stderr``.
  :param k0: Adjust the prefactor for the Metropolis rule when calculating rates.
  :param temp: Specify temperature in Celsius.
  :type name: string
  :type seq: string
  :type sfile: string
  :type temp: float
  :type k0: float
  :type force: bool
  :type verb: bool

  :return: A list of produced files containing ``barriers`` results. \
      [sfile, bfile, rfile, efile, pfile]
  :rtype: list
  """

  if which(barriers) is None :
    print barriers, "is not executable"
    print """ 
    You need to install *barriers*, which you can download from the 
    ViennaRNA package homepage: http://www.tbi.univie.ac.at/RNA/Barriers/
    
    If you have installed the program, make sure that the path you specified 
    is executable.
    """
    raise RuntimeError('Could not find executable')

  if not sfile or not os.path.exists(sfile) : 
    sfile = sys_suboptimals(name, seq, 
        verb=verb, noLP=noLP, circ=circ, temp=temp, force=force)

  bfile = name+'.bar'
  efile = name+'.err'
  rfile = name+'.rts'
  psfile= name+'.ps'

  if not force and \
      os.path.exists(bfile) and \
      os.path.exists(rfile) and \
      os.path.exists(psfile) and \
      os.path.exists(efile):
      if verb : print >> sys.stderr, bfile, "<= Files exist"
      return [sfile, bfile, efile, rfile, psfile]

  if gzip :
    barcall = ['zcat', sfile, '|', barriers]
  else :
    barcall = ['cat', sfile, '|', barriers]

  if noLP : 
    barcall.extend(['-G', 'RNA-noLP']) 
  else :
    barcall.extend(['-G', 'RNA'])

  if moves == 'single-base-pair' :
    barcall.extend(['--moves=noShift'])
  elif moves == 'shift' :
    barcall.extend(['--moves=Shift'])
  else :
    print >> sys.stderr, 'Invalid move-set in sys_barriers()'

  # buggy barriers
  if subopt_reaches_minh(sfile, minh) : 
    barcall.extend(["--minh", str(minh)])
    barcall.extend(["--max", str(int(maxn))])
  barcall.extend(["-T", str(temp)])

  if rates : barcall.extend(['--rates']) 
  if bsize : barcall.extend(['--bsize']) 
  if saddle: barcall.extend(['--saddle']) 

  if mfile : barcall.extend(["--mapstruc", mfile])

  if verb :
    print >> sys.stderr, ' '.join(barcall), '2>', efile, '>', bfile

  with open(sfile, 'r') as shandle, \
      open(bfile, 'w') as bhandle, \
      open(efile, 'w') as ehandle:
    process = sub.Popen(' '.join(barcall),
        stdin=shandle,stdout=bhandle,stderr=ehandle, shell=True)
    process.communicate(None)
    if process.returncode :
      print >> sys.stderr, \
          "ERROR: process terminated with return code", process.returncode
      print >> sys.stderr, 'Error:', ' '.join(barcall), '2>', efile, '>', bfile
      raise RuntimeError

  if rates :
    if k0 != 1.0:
      with open('rates.out') as prerates, open(rfile, 'w') as rates:
        for line in prerates:
          newline = "".join(map(
            "{:10.4g}".format, [float(x)*k0 for x in line.strip().split()]))
          rates.write(newline+"\n")
      os.remove('rates.out')
    else :
      os.rename('rates.out', rfile)
    os.remove('rates.bin')
    os.remove('treeR.ps')
  os.rename('tree.ps', psfile)

  return [sfile, bfile, efile, rfile, psfile]

def sys_suboptimals(name, seq, 
    RNAsubopt='RNAsubopt',
    ener=None,
    temp=37.0,
    verb=False,
    noLP=False,
    circ=False,
    opts=[],
    sort=['|', 'sort', '-T', '/tmp', '-k3r', '-k2n'],
    gzip=True,
    force=False):
  """ **Perform a system-call of the program ``RNAsubopt``.**

  The print the results into a file and return the filename. This wrapper will
  produce two output files from ``STDIN`` and ``STDERR``, respectively. 

  :param name: Name of the sequence used to name the output file.
  :param seq: Nucleic acid sequence.
  :param RNAsubopt: path to executable
  :param ener: Specify energy range
  :param temp: Specify temperature in Celsius.
  :param verb: Print the current system-call to ``stderr``.
  :param noLP: exclude lonely-base-pairs in suboptimal structures
  :param circ: compute density of states
  :param opts: more options for ``RNAsubopt``

  :param force: Overwrite existing files with the same name.

  :type name: string
  :type seq: string
  :type RNAsubopt: string
  :type ener: float
  :type temp: float
  :type noLP: bool
  :type circ: bool
  :type opts: list
  :type force: bool
  :type verb: bool

  :return: Fliename of the file containing ``RNAsubopt`` results
  :rtype: string
  """

  if which(RNAsubopt) is None :
    print RNAsubopt, "is not executable"
    print """ 
    You need to install *RNAsubopt*, which is part of the ViennaRNA
    package, download: http://www.tbi.univie.ac.at/RNA 
    
    If you have installed the program, make sure that the path you specified 
    is executable.
    """
    raise RuntimeError('Could not find executable')

  if ener is None :
    ener, nos = sys_subopt_range(seq, verb=verb, 
        RNAsubopt=RNAsubopt, noLP=noLP, circ=circ, temp=temp)
    if verb: print >> sys.stderr, \
        "Energy-Update: {:.2f} kcal/mol to compute {} sequences".format(ener, nos)

  sfile = name+'.spt.gz' if gzip else name + '.spt'

  if os.path.exists(sfile) and not force : 
    if verb: print >> sys.stderr, sfile, "<= File exists"
    return sfile

  sptcall = [RNAsubopt, "-e {:.2f} -T {:.2f}".format(ener, temp)]
  if noLP : sptcall.append("--noLP")
  if circ : sptcall.append("--circ")
  sptcall.extend(opts)
  sptcall.extend(sort)

  if gzip : sptcall.extend(('|', 'gzip', '--best')) 

  if verb :
    print >> sys.stderr, 'echo', '"'+seq+'"', '|', ' '.join(sptcall), '>', sfile

  with open(sfile, 'w') as shandle:
    Psubopt = sub.Popen([' '.join(sptcall)],
        stdin=sub.PIPE, stdout=shandle, shell=True)
    Psubopt.communicate(seq)
    if Psubopt.returncode :
      print >> sys.stderr, \
        "ERROR: process terminated with return code", Psubopt.returncode
      print >> sys.stderr, 'Error:', \
          'echo', '"'+seq+'"', '|', ' '.join(sptcall), '>', sfile
      raise RuntimeError
  return sfile

def sys_subopt_range(seq,
    RNAsubopt='RNAsubopt',
    nos=5100000,
    maxe=30.0,
    temp=37.0,
    noLP=False,
    circ=False,
    verb=False):
  """ Compute an energy range that computes a given number of structures.
  
  .. note:: If your RAM is big enough, ``barriers`` can be compiled to read \
    billions of structures, however computations with more than 10.000.000 \
    structures may take a looong time.

  :param seq: nucleic acid sequence
  :param RNAsubopt: path to executable
  :param nos: number of structures
  :param maxe: an upper bound of the energy range
  :param temp: temperature in Celsius
  :param noLP: exclude lonely-base-pairs in suboptimal structures
  :param circ: compute density of states
  :param verb: print verbose output to ``stderr``
  
  :return: (energy-range, number-of-structures)
  :rtype: tuple
  """

  if which(RNAsubopt) is None :
    print RNAsubopt, "is not executable"
    print """ 
    You need to install *RNAsubopt*, which is part of the ViennaRNA
    package, download: http://www.tbi.univie.ac.at/RNA 
    
    If you have installed the program, make sure that the path you specified 
    is executable.
    """
    raise RuntimeError('Could not find executable')


  num, nump = 0, 0
  e = maxe if nos is 0 else 2.
  ep = e-1
  while (num < nos+1) :
    if verb: 
      print "Energy: ", "{:.2f}".format(float(e))

    sptcall = [RNAsubopt, "-D -e {:.2f} -T {:.2f}".format(float(e), temp)]
    if circ : sptcall.append("--circ")
    if noLP : sptcall.append("--noLP")

    process = sub.Popen([' '.join(sptcall)], 
            stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
    output, err = process.communicate(seq)
    if err : print "ERROR:", err

    structures = 0
    for l in output.split("\n")[1:-1]:
      [interval,value] = l.split()
      structures += int(value)

      # Make sure that nump fits to the preset ep value
      if nump == 0 :
        if ep*10 == int(interval): 
          nump = structures 

      # Break the loop if you reach the number of structures in the output
      if (nos and structures >= nos - (nos * 0.01)) or float(interval)/10 >= maxe:
        #print "break", maxe, nos, structures, "E:", float(interval)/10
        e   = float(interval)/10
        num = structures
        break # return to while

    else : # if the last loop didn't break, duh!
      num=structures
      if num > nos or num == nump : 
        e = float(interval)/10
        break

      new = e + (e-ep)/(math.log(float(num)/nump))*(math.log(float(nos)/num))
      ep, nump = e, num
      e = new if new < maxe else maxe
      if abs(e-ep)<0.1 : 
        e=ep
        break
      continue
    break # end while after for loop exit

  return e, num

# ***************** #
# private functions #
# ................. #

def subopt_reaches_minh(fname, minh):
  """ Internal function to report on whether the energy-range of suboptimal
  structures exceeds the ``--minh`` options for ``barriers``. If this is not
  the case, the ``--minh`` option can cause segmentation faults.

  :param fname: Filename of **gzipped and sorted** RNAsubopt result 
  :param minh: The ``barriers --minh`` value

  :type fname: string
  :type minh: float

  :return: True or False
  """
  with gzip.open(fname, 'r') as f:
    for i, l in enumerate(f):
      if i == 0:
        continue
      elif i == 1:
        mfe = l.strip().split()[1] 
      else:
        sen = l.strip().split()[1] 
        if float(sen) - float(mfe) > minh :
          return True
  return False

def which(program):
  """ Emulates the unix ``which`` command. Snatched from:
    `http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python`

    :param program: executable
    :type program: string

    :returns: path-to-executable or None
  """
  def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

  fpath, fname = os.path.split(program)
  if fpath:
    if is_exe(program):
      return program
  else:
    for path in os.environ["PATH"].split(os.pathsep):
      path = path.strip('"')
      exe_file = os.path.join(path, program)
      if is_exe(exe_file):
        return exe_file

  return None

def main():
  import ribolands.utils as rnu
  
  # Quick set test model params """
  name, seq = rnu.parse_vienna_stdin()
  print name, seq

  #param='RNA'
  #dangle='some'
  s_nos = 1e6
  s_maxe = 20
  s_ener = None
  s_patch=[] #['|', 'pylands_spatch.py', '--theo']

  b_minh=2.0
  b_maxn=20

  temp=37.0
  circ=False
  noLP=True
  verb=True
  force=True

  if s_ener is None :
    s_ener, s_nos = sys_subopt_range(seq, nos=s_nos, maxe=s_maxe, verb=verb)
  sfile = sys_suboptimals(name, seq, 
      ener=s_ener, 
      noLP=noLP,
      opts=s_patch,
      verb=verb, 
      force=force)

  [sfile, bfile, efile, rfile, psfile] = sys_barriers(name, seq, sfile, 
      minh=b_minh, maxn=b_maxn, rates=True, verb=verb, noLP=noLP, force=force)
  tfile = sys_treekin(name, seq, bfile, rfile, 
      p0=['2=1'], t0=1e-6, ti=1.02, t8=1e10, verb=verb, force=force)
  pfile = rnu.plot_simulation(name, seq, tfile, 
      ylim=(0,1), xlim=(1e-2, 1e10), lines=[], force=force)

  # BCG = barriersCG(mfile, efile)
  RM = rnu.parse_ratefile(rfile)
  BT = rnu.parse_barfile(bfile, seq=seq)

  print RM, BT
  print sfile, bfile, efile, rfile, psfile, tfile, pfile

  return

if __name__ == '__main__':
  main()

