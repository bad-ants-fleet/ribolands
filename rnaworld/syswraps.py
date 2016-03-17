#!/usr/bin/env python
#
# Written by Stefan Badelt <stef@tbi.univie.ac.at>
#
# systemcalls of RNAsubopt, barriers and treekin
# tested only under linux
#

import os
import re
import sys
import gzip
import math
import subprocess as sub

""" start of public functions """
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
  """
  Do the treekin simulation and print the results after each simulation as
  finished to stdout.

  :return: None
  """
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
        sys.exit('over and out')

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
  """ single barriers run, produce missing files """

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
      raise SystemExit

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
  """ Call RNAsubopt """

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
      raise SystemExit
  return sfile

def sys_subopt_range(seq,
    RNAsubopt='RNAsubopt',
    temp=37.0,
    verb=False,
    noLP=False,
    circ=False,
    maxe=30.0,
    nos=5100000):
  """ Compute energy range for given number of structures """

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

""" start of internal functions """
def subopt_reaches_minh(fname, minh):
  """ buggy barriers """
  with gzip.open(fname, 'r') as f:
    for i, l in enumerate(f):
      if i == 0:
        continue
      elif i == 1:
        mfe = l.strip().split()[1] 
      else:
        sen = l.strip().split()[1] 
        if float(sen) - float(mfe) > minh :
          return 1
  return 0


