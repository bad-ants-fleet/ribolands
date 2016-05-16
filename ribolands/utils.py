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

import re
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

