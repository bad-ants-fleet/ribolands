#!/usr/bin/env python

__version__ = "0.1"
__author__  = "Stefan Badelt"
__email__   = "stef@tbi.univie.ac.at"

"""
  Style: If you use strings as comments, like here, use double quotes, such
  that single quotes are available for uncommenting large parts during testing

  *) do not exceed 80 characters per line
  *) indents: 2x whitespace, no tab characters!

  vim-config settings:
  set textwidth=80
  set tabstops ...
"""

import re
import sys
import os
import argparse
import collections as c
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

""" * * * * * * TODO * * * * * *
  *) make a python module out of this 

      ViennaRNA
      -> import RNA as vrna

      Python Libraries
      -> vrna.landscapes
      -> vrna.landscapes.syscalls
      -> vrna.landscapes.walks

      -> vrna.landscapes.development

      Scripts/Tests:
      -> BarMap.py
      -> DrTransformer.py
      -> spatch.py
      -> interkin.py

  *) on import:
        -> catch errors if RNA, RNAsubopt, barriers and treekin don't exist,
          warn otherwise
  *) make barriers work in its own temporary directory 
    and *then* rename the standard file names (rates.out)
    alternative: warn that the user must not make parallel computations

  *) use function wrappers (?) to specify a common set of model details
  *) think of naming I/O files (*.bar, *.spt, *.spt.gz, *.tkn, *.rts, *.ps, *.pdf ...)
  *) compare spatch-hack with interkin 
  *) write test files
  *) sphinx documentation
"""

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

'''
def cofold_barriers(_name, species,
    spatch=['|', '~/WORK/Python/spatch.py'],
    ener = None,
    bionly=False, #DEBUGGING
    barriers='barriers',
    minh=0.001,
    maxn=50,
    k0=1.0,
    temp=37.0,
    moves='single-base-pair',
    gzip=True,
    rates=True,
    bsize=False,
    saddle=False,
    force=False,
    verb=False):
  """ Simulate inter-molecular folding
  Compute folding kinetics 
    returns the biggest connected component and a list of nodes
    TODO: 
      *) @homo-dimers: RNAsubopt-barriers does not correct for two-fold symmetry!
  """
  
  import networkx as nx

  noLP=False # MUST BE FALSE!
  #minh=0.5  # MUST BE < 1 (?)
  circ=False # MUST BE FALSE!
  mfile=''

  RG = nx.MultiDiGraph()
  sortednodes = c.OrderedDict()
  hypercount  = [0]

  for spe, seq in species :
    name = _name + "_" + spe

    if '&' in name :
      name = name.translate(string.maketrans("&", "_"))

      """ additional true dimer only run """
      if not bionly :
        spatch.append('-d')

        if ener is None : sys.exit('specify energy range!')
        sfile = sys_suboptimals(name, seq, 
            ener=ener, opts=spatch, verb=verb, force=force)
        [sfile, bfile, efile, rfile, psfile] = sys_barriers(name, seq, sfile, 
            barriers=barriers, minh=minh, maxn=maxn, k0=k0, temp=temp, 
            verb=verb, force=force)

        # this function modifies sortednodes and hypercount... not nice!
        add_edges(RG, bfile, rfile, 'uni', spe, seq, sortednodes, hypercount)

        args.s_patch.pop()

      [bfile, efile, rfile, ps] = bm.barriers(seq, args)
      if bionly :
        add_edges(RG, bfile, rfile, 'bionly', spe, seq, sortednodes, hypercount)
      else :
        add_edges(RG, bfile, rfile, 'bi', spe, seq, sortednodes, hypercount)
    else :
      if not bionly :
        [bfile, efile, rfile, ps] = bm.barriers(seq, args)
        add_edges(RG, bfile, rfile, 'uni', spe, seq, sortednodes, hypercount)
    args.name = name

    print "connected components", 
    print [len(comp) for comp in 
        sorted(nx.strongly_connected_components(RG), key=len, reverse=True)]

  for e, comp in enumerate(
      sorted(nx.strongly_connected_components(RG), key=len, reverse=True)) :
    #print comp
    if e == 0 : continue
    print "Discarding connected component!", len(comp), comp
    for i in comp: sortednodes[i] = 0

  return max(nx.strongly_connected_component_subgraphs(RG), key=len), sortednodes


def add_edges(RG, bfile, rfile, mode, spe, seq, snodes, hcount):
  """ Add the transition edges, no symmetry correction! """
  RM = []
  lmins = [] #map(lambda(x): (spe, x[1]), nal.parse_barfile(bfile, seq=seq))
  for line in nal.parse_barfile(bfile, seq=seq):
    lmins.append((spe, line[1]))
    if (spe, line[1]) not in snodes :
      snodes[(spe,line[1])] = 0
  with open(rfile) as rates :
    for line in rates :
      """ Can we correct rates for two-fold rotational symmetry ? """
      RM.append((map(float, line.strip().split())))
  pDF = pd.DataFrame(RM, columns=lmins, index=lmins)

  if mode == 'uni' :
    for col in pDF : 
      for row in pDF :
        if row == col : continue
        if pDF[col][row] == 0. : continue
        if pDF[row][col] == 0. :
          print "backward rate missing!!!"
        #The forward rate row => col is in RM[col][row]
        #print "Adding:", row, col, pDF[col][row]
        RG.add_weighted_edges_from([(row, col, pDF[col][row])])
        snodes[row]=1
        snodes[col]=1
  elif mode == 'bi' or mode == 'bionly' :
    #fakemin
    for col in pDF : 
      for row in pDF :
        if row == col : continue
        if pDF[col][row] == 0. : continue
        #The forward rate row => col is in RM[col][row]
        cplR = get_complexes(row)
        cplC = get_complexes(col)
        # good enough to identify pairwise interactions:
        if len(cplR) == len(cplC) : 
          """ WARNING: adding unimolecular rate from cofold tree """
          if mode == 'bi' : continue
          if len(cplR) == 1 :
            """ Dimer to Dimer """
            RG.add_weighted_edges_from([(row, col, pDF[col][row])])
            snodes[row]=1
            snodes[col]=1
          elif len(cplR) == 2 :
            """ 2xMon to 2xMon """
            if cplR[0] != cplC[0] and cplR[1] != cplC[1] : 
              """ this fucks up simulations, skip it! """
              continue
            if cplR[0] != cplC[0] :
              #print "adding", cplR[0], cplC[0]
              edict = RG.get_edge_data(cplR[0], cplC[0])
              if edict :
                if edict[0]['weight'] < pDF[col][row] :
                  edict[0]['weight'] = pDF[col][row]
              else :
                RG.add_weighted_edges_from([(cplR[0], cplC[0], pDF[col][row])])
              snodes[cplR[0]]=1
              snodes[cplC[0]]=1
            if cplR[1] != cplC[1] :
              #print "adding", cplR[1], cplC[1], pDF[col][row], row, col 
              edict = RG.get_edge_data(cplR[1], cplC[1])
              if edict :
                if edict[0]['weight'] < pDF[col][row] :
                  edict[0]['weight'] = pDF[col][row]
              else :
                RG.add_weighted_edges_from([(cplR[1], cplC[1], pDF[col][row])])
              snodes[cplR[1]]=1
              snodes[cplC[1]]=1
          else :
            sys.exit('no support for compexes of size > 2')
        else :
          hcount[0] += 1
          hnode = ('HYPER', str(hcount[0]))
          for tup in cplR+cplC :
            if not RG.has_node(tup):
              #print tup, "recovered!"
              snodes[tup]=1
          #print "Adding:", row, col, pDF[col][row]
          for tup in cplR :
            RG.add_weighted_edges_from([(tup, hnode, pDF[col][row])])
          for tup in cplC :
            RG.add_weighted_edges_from([(hnode, tup, pDF[col][row])])
  else :
    print """ specify valid modus (uni/bi) """
  
  return
'''

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
  spatch=['|', '~/WORK/Python/spatch.py', '--theo']

  verb=True

  #ener, nos = sys_subopt_range(seq, nos=1000000, maxe=20, verb=verb)
  sfile = sys_suboptimals(name, seq, 
      ener=18, 
      noLP=noLP,
      opts=spatch,
      verb=verb, 
      force=force)

  [sfile, bfile, efile, rfile, psfile] = sys_barriers(name, seq, sfile, 
      minh=2.0, maxn=20, rates=True, verb=verb, noLP=noLP, force=force)
  tfile = sys_treekin(name, seq, bfile, rfile, 
      p0=['2=1'], t0=1e-6, ti=1.02, t8=1e10, verb=verb, force=force)
  pfile = plot_simulation(name, seq, tfile, 
      ylim=(0,1), xlim=(1e-2, 1e10), lines=[], force=force)

  # BCG = barriersCG(mfile, efile)
  RM = parse_ratefile(rfile)
  BT = parse_barfile(bfile, seq=seq)

  print RM, BT
  print sfile, bfile, efile, rfile, psfile, tfile, pfile

  return

if __name__ == '__main__':
  main()

