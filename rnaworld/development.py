#!/usr/bin/env python

"""
  Coded by: Stefan Badelt <stef@tbi.univie.ac.at>
  University of Vienna, Department of Theoretical Chemistry

  -*- Style -*- 
  Use double quotes or '#' for comments, such that single quotes are available
  for uncommenting large parts during testing

  *) do not exceed 80 characters per line
  *) indents: 2x whitespace, no tab characters!

  -*- VIM config -*- 
  set textwidth=80
  set ts=2 et sw=2 sts=2
"""

import re
import sys
import os
import argparse
import collections as c
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable



""" Thinktank for a new and better interface:

  ******** Example *********
  
  context = pyrna_environment(tempC=37.0, NaCL=1.0, theo=False, ...)
  model   = pyrna_energy_model(params='turner', dangle=2, noGU=False, tetraloop=True)
  
  romeo = pyrna_nucleic_acid(seq, circ=False, environment=context, emodel=model)
  
  romeo.set_environment.tempC(37.0)
  romeo.set_environment.kT(0.61)
  romeo.set_energy_model.params('turner1999')
  
  [ss, en] = romeo.mfe(tempC=15, circ=True) # Sets the variable for this run only
  # alternative: pyrna_mfe(romeo, tempC=15, circ=True)
  efe = romeo.efe
  pfc = romeo.pfc
  ffile = romeo.sys_fold
  
  romeo.energy_model.noLP(True)
  romeo.energy_model.k0(2e5)
  sfile = romeo.sys_suboptimals(fname=romeo, ener=None, maxn=5000, maxe=30.00)
  [bfile, rfile, ... ] = romeo.sys_barriers(minh=1, maxn=50, rates=True, k0=1)
  tfile = romeo.sys_treekin(t0=1e-6, ti=1.02, t8=1e6)
  
  what about two sequences?
  context = romeo.get_environment()
  model = romeo.get_energy_model()
  julia = pyrna_nucleic_acid(seq, context=context, model=model)
  
  
  couple = romeo + julia # only works if context is the same (?)
  
  cut_point = couple.get_cut_points
  couple.mfe(cut_point=-1)
  copule.sys_suboptimals
  
  onemore += felix # will be prohibited, unless using NUPACK

  ******** Structure *********

  Initialize parent objects:
    * vrna_environment # temp, circ
    * energy_model # 
    * model_options(subopt, barriers, treekin)
  
  Initialize RNA object and inherit from parents, if they are set:
    * pyrna_nucleic_acid(seq, env=env, emod=emod)
  
  Call functions with obligatory one-time parameters

"""

""" under construction """
class pyrna_context:
  """ Set model details to standard vrna details """
  temperature = 37.0
  theophylline = None
  tetracycline = None
  ions = None

# inhertiance: children have access to parent methods
class pyrna_molecule: 
  """ Define basic properties of the molecule. 
    This is the partent class
  
  """
  def __init__(self, sequence):
    self.sequence = sequence
    circ = False
    print "Creating an RNA object"

  def energy_model(self):
    parameter_file = 'turner1999'
    dangles = 2 # [0: none, 1: ?, 2: some, 3: mathews]
    GUclosing = True
    tetraloopen = True


class pyrna_options:
  def __init__(self):
    pass

  def general(self):
    noLP = False

  def suboptimals(self):
    subopt_energy = None
    maxspt_energy = 30.0
    subopt_structs= 5000

  def barriers(self):
    rates=True
    bsize=False
    saddle=False
    minh = 0.01
    maxn = 100

  def kinetics(self):
    t0 = 1e-6
    ti = 1.02
    t8 = 1e10


""" Under Testing """

class cd:
  """ Context manager for changing the current working directory 
    http://stackoverflow.com/questions/431684/how-do-i-cd-in-python 
  """
  def __init__(self, newPath):
    self.newPath = os.path.expanduser(newPath)

  def __enter__(self):
    self.savedPath = os.getcwd()
    os.chdir(self.newPath)

  def __exit__(self, etype, value, traceback):
    os.chdir(self.savedPath)

'''
with cd("~/Library"):
  s.call("ls")
'''
def model_details(force=False, tmp=27, **kwars):
  def wrapper(name):
    return "bla-{0}".format(func(name))
  return wrapper


""" Main public functions """

def plot_simulation(name, seq, tfile, 
    title='',
    plim=1e-2,
    lines=[],
    ylim=(None,None),
    xlim=(None,None), 
    verb=True, 
    lilog=None,
    force=True):
  """ Plot a (particular) list of trajectories
  t [x x x x x]
  1 [y y y y y]
  2 [y y y y y]
  3 [y y y y y]
  4 [y y y y y]
  """
  lines=set(lines)
  title = title if title else name

  nxy=[]
  with open(tfile) as tkn :
    for line in tkn :
      if re.match('#', line): 
        continue
      nxy.append(map(float, line.strip().split()))

  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  if None not in ylim : ax.set_ylim(ylim)
  if None not in xlim : ax.set_xlim(xlim) 

  if lilog :
    ''' Make the second part of the plot logarithmic'''
    ax.set_xscale('linear')
    ax.set_xlim((0, lilog))
    divider = make_axes_locatable(ax)
    axLog = divider.append_axes("right", size=2, pad=0.0, sharey=ax)
    axLog.set_xlim((lilog, 1e5))
    axLog.set_xscale('log')
    #axLog.set_yticklabels([''])
  else :
    ax.set_xscale('log')

  time = []
  for e, traject in enumerate(zip(*nxy)):
    if e==0: 
      time = traject
      continue
    if lines and e not in lines : continue
    if plim and max(traject) < plim : continue
    if ylim and max(traject) > ylim : continue
    p, = ax.plot(time, traject, '-')
    if lilog: p, = axLog.plot(time, traject, '-')
    p.set_label("lmin {:d}".format(e))
    #p.set_color(axcolors[finalmin])

  fig.set_size_inches(7,3)
  fig.text(0.5,0.95, title, ha='center', va='center')
  plt.xlabel('time [seconds]', fontsize=11)
  plt.ylabel('occupancy [mol/l]', fontsize=11)

  # """ Add ticks for 1 minute, 1 hour, 1 day, 1 year """
  # plt.axvline(x=60, linewidth=1, color='black', linestyle='--')
  # plt.axvline(x=3600, linewidth=1, color='black', linestyle='--')
  # plt.axvline(x=86400, linewidth=1, color='black', linestyle='--')
  # plt.axvline(x=31536000, linewidth=1, color='black', linestyle='--')
  plt.legend()

  pfile = name+'.pdf'
  plt.savefig(pfile, bbox_inches='tight')

  return pfile

""" utilities """


def barriersCG(mfile, efile, verb=False):
  """ Takes a sequence and a file with structures, populations, ...
  (mfile). Computes a barrier tree for the given sequence and returns 
  a dictionary with the correctly mapped populations:
  (requires barriers-v1.6 or later)

  mfile in form:
  structure energy population (+ lost columns)
  ...
  
  :return: Dictionary[lmin] = population
  """
  
  BarCG = c.defaultdict(float)
  with open(mfile, 'r') as m, open(efile, 'r') as mapstruc:
    for line in mapstruc:
      if re.match('[\.\(\)]+', line.strip()):
        gstr, sptidx, energy, fmin, fminT, gmin, gminT = line.strip().split()
        [old, pop, en] = m.readline().strip().split()[0:3]
        if verb : print [gminT, gstr, old, en, pop]
        BarCG[int(gminT)] += float(pop)
      elif re.match('not in hash', line.strip()):
        [old, pop, en] = m.readline().strip().split()[0:3]
        print >> sys.stderr, old, en, "structure not in hash"
        sys.exit('over and out')
      elif re.match('not yet assigned', line.strip()):
        [old, pop, en] = m.readline().strip().split()[0:3]
        print >> sys.stderr, old, en, "structure not yet assigned"
        sys.exit('over and out')
  return BarCG


""" A little example for main """

def main():
  import rnaworld.syswraps as rns
  import rnaworld.utils as rnu

  """ xxx """
  parser = argparse.ArgumentParser()
  parser.add_argument("-p","--path", 
    help="Specify path to local file or directory", 
    default='/dev/null',
    action = 'store')
  parser.add_argument("-n","--name",
    help="Name your output", 
    default= '')
  parser.add_argument("-i","--iterations",
    help="Do stuff multiple times",
    type=int,
    default=5) 
  parser.add_argument("-v","--verbose",
    help="Verbose output",
    action="store_true")
  args = parser.parse_args()
  
  """ Model Parameters """
  name, seq = read_vienna_stdin()
  print name, seq
  param='RNA'
  dangle='some'
  temp=37.0
  circ=False
  noLP=True
  force=True
  spatch=['|', 'pylands_spatch.py', '--theo']

  verb=True

  #ener, nos = sys_subopt_range(seq, nos=1000000, maxe=20, verb=verb)
  sfile = rns.sys_suboptimals(name, seq, 
      ener=18, 
      noLP=noLP,
      opts=spatch,
      verb=verb, 
      force=force)

  [sfile, bfile, efile, rfile, psfile] = rns.sys_barriers(name, seq, sfile, 
      minh=2.0, maxn=20, rates=True, verb=verb, noLP=noLP, force=force)
  tfile = rns.sys_treekin(name, seq, bfile, rfile, 
      p0=['2=1'], t0=1e-6, ti=1.02, t8=1e10, verb=verb, force=force)
  pfile = plot_simulation(name, seq, tfile, 
      ylim=(0,1), xlim=(1e-2, 1e10), lines=[], force=force)

  # BCG = barriersCG(mfile, efile)
  RM = rnu.parse_ratefile(rfile)
  BT = rnu.parse_barfile(bfile, seq=seq)

  print RM, BT
  print sfile, bfile, efile, rfile, psfile, tfile, pfile

  return

if __name__ == '__main__':
  main()

