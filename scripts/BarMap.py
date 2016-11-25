#!/usr/bin/env python

""" 
  BarMap 2 -- cotranscriptional folding with *barriers* and *treekin*

  :requires: python-v.2.7, RNAsubopt, barriers-v1.6, treekin-v0.4

  use two whitespace characters as tab when editing this file!
  vim-config = set: ts=2 et sw=2 sts=2
"""

import os
import re
import sys
import argparse
import numpy as np
import subprocess as s # barmap_treekin()
import collections as c
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import ribolands as ril
from ribolands.syswraps import which

def barmap_plot(name, seq, tfiles, plist, args):
  """ Description, 
    TODO 
  """
  verb = args.verbose
  start= args.start
  stop = args.stop
  title= args.plot_title if args.plot_title else name
  plim = args.plot_cutoff

  t8   = args.t8
  tX   = args.tX
  linlog= args.plot_linlog

  lin_time = 0.0
  tot_time = 0.0

  reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')
  # http://www.regular-expressions.info/floatingpoint.html

  tlist = [[] for i in range(len(plist)+1)]
  for e, l in enumerate(range(start,stop+1)):
    tfile = tfiles[e]

    ''' Get raw data => nxy
     t 1 2 3 4 
    [x y y y y]
    [x y y y y] '''
    nxy=[]
    with open(tfile) as tkn :
      for line in tkn :
        if re.match('#', line): 
          continue
        elif not reg_flt.match(line):
          nxy.append([t8, 1.0])
          if verb :
            print >> sys.stderr, "WARNING:", tfile, "does not contain data!"
          break
        nxy.append(map(float, line.strip().split()))

    ''' Make list of trajectories
    t [x x x x x]
    1 [y y y y y]
    2 [y y y y y]
    3 [y y y y y]
    4 [y y y y y]

    and use the pathlist with lmin information:
    [0 0 0 3 3 2 2 2 0]
    [1 1 1 5 5 6 8 9 0]
    [3 3 3 2 2 1 1 2 1]
    [2 2 2 2 2 1 1 2 1]
    [0 9 9 9 9 9 1 2 1]
    ...
    to replace every lmin at the current transcription step with the
    trajectory (or None) => avoid duplicates with seen=set()!
    '''
    datalength = len(np.array(nxy)[:,0])
    timeline = np.array(nxy)[:,0] + tot_time
    lin_time = tot_time
    tot_time = timeline[-1]
    tlist[0].extend(timeline)
    seen = set()
    for e, idx in enumerate(np.array(plist)[:,l-start]):
      # Make sure that we are above the limit of detection
      if idx != 0 and idx not in seen and max(np.array(nxy)[:,idx]) > plim :
      #if idx != 0 and max(np.array(nxy)[:,idx]) > plim:
        line = np.array(nxy)[:,idx]
      else :
        line = [None for i in range(datalength)]
        seen.add(idx)
      tlist[e+1].extend(line)

  lin_time = (stop-start)*float(t8)
  if verb :
    print >> sys.stderr, "Total simulation time:", tot_time
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  #ax.set_ylim([0,1.01])
  ax.set_xscale('linear')

  ''' Make the second part of the plot logarithmic'''
  if linlog :
    ax.set_xlim((0, lin_time))
    divider = make_axes_locatable(ax)
    axLog = divider.append_axes("right", size=2.5, pad=0.0, sharey=ax)
    axLog.set_ylim([0,1.01])
    #axLog.set_yticklabels([''])
    axLog.set_xlim((lin_time, tX))
    axLog.set_xscale('log')

  axcolors = c.defaultdict(int)
  for e, traject in enumerate(tlist):
    if e==0: continue
    if max(tlist[e]) < plim : continue
    finalmin = plist[e-1][-1]
    '''color='b' will be overwritten, but somehow its necessary'''
    lin, = ax.plot(tlist[0], tlist[e], '-', color='b')
    if linlog :
      log, = axLog.plot(tlist[0], tlist[e])
    if axcolors[finalmin] == 0 :
      lin.set_label("lmin {:d}".format(finalmin))
      if linlog :
        log.set_label("lmin {:d}".format(finalmin))
      axcolors[finalmin] = ax._get_lines.color_cycle.next()
    lin.set_color(axcolors[finalmin])
    if linlog :
      log.set_color(axcolors[finalmin])
   
  plt.legend()
  fig.text(0.5,0.95, title, ha='center', va='center')
  fig.text(0.5,0.04, 'time [seconds]', ha='center', va='center')
  fig.text(0.06,0.5, 'occupancy', ha='center', va='center', rotation='vertical')
  pfile = "{}.pdf".format(name)
  plt.savefig(pfile)

  return pfile

def barmap_treekin(bname, seq, bfiles, plist, args):
  """ write treekin files from bfiles and plist information """
  tmpdir= args.tmpdir
  start = args.start
  stop  = args.stop 
  verb  = args.verbose
  tX    = args.tX
  t8    = args.t8
  cutoff= args.cutoff

  p0 = args.p0
  tt = 0
  tfiles = []
  reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')
  # http://www.regular-expressions.info/floatingpoint.html

  for e, l in enumerate(range(start,stop+1)):
    t8 = tX if l == stop else t8

    cseq  = seq[0:l]
    cname = "{}-t8_{}-len_{}".format(bname,t8,l)
    [bfile, efile, rfile, psfile] = bfiles[e]

    ctfile, _ = ril.sys_treekin(cname, cseq, bfile, rfile, 
        treekin = args.treekin,
        repl=None,
        #useplusI=True,
        p0=p0, 
        t0=args.t0, 
        ti=args.ti, 
        t8=t8, 
        verb=verb, 
        force=args.force)
    tfiles.append(ctfile)

    lastlines = s.check_output(['tail', '-2', ctfile]).strip().split("\n")
    if reg_flt.match(lastlines[0]):
      tt += float(lastlines[0].split()[0])
      if l < stop :
        curlmin = np.array(plist)[:,e]
        newlmin = np.array(plist)[:,e+1]
        p0 = set_p0(bfile, l, lastlines, curlmin, newlmin,
            cutoff, verb)
      else :
        for (i, pop) in enumerate(lastlines[0].split()):
          if i != 0 and float(pop) > cutoff :
            print "{:3d} {:3d} {:f} {:s}".format(l, i, 
                float(pop), get_structure(bfile, i))
    else :
      tt += t8
      print "{:3d} {:3d} {:f} {:s}".format(l, 1, 1.0, get_structure(bfile, 1))
      if get_structure(bfile, 2):
        print >> sys.stderr, \
            "No valid time-course in {:s}: {:s}".format(ctfile, lastlines[0])
        with open(efile, 'r') as err:
          print >> sys.stderr, err.read().strip()
        sys.exit('over and out')

  return tfiles

def barmap_mapping(_bname, seq, args):
  """ Parse mapping info into pathlist """
  tmpdir= args.tmpdir
  start = args.start
  stop  = args.stop 
  verb  = args.verbose
  force = args.force

  plist = []
  mapinfo = "{}-map_{:d}_to_{:d}.map".format(_bname, start, stop)

  if os.path.exists(mapinfo) and not force:
    if verb :
      print >> sys.stderr, mapinfo, "<= File exists"
    with open(mapinfo, 'r') as m:
      m.readline().strip()
      for line in m:
        path = [int(x) for x in line.strip().split(" => ")]
        plist.append(path)
  else :
    mlist = []
    for l in range(start, stop+1):
      cseq  = seq[0:l]
      cname = "{}-len_{:02d}".format(_bname,l)
      pname = "{}-len_{:02d}".format(_bname,l-1)

      if os.path.exists(pname + '.bar'):
        if verb : 
          print >> sys.stderr, "Get mapping info {:d} -> {:d}".format(l-1,l)
        mlist.append(get_mapping_dict(pname + '.bar', cname + '.err'))
 
    plist = pathlist(mlist)
    with open(mapinfo, 'w') as m:
      m.write("# Mapping information for barrier trees {:d} to {:d}\n".format(
        start, stop))
      for path in plist: 
        m.write(" => ".join(map("{:4d}".format, path)) + "\n")
  return plist

def barmap_barriers(_bname, seq, sfiles, args):
  """ Print barriers files and mapping info """

  bfiles = []
  for e, l in enumerate(range(args.start,args.stop+1)):
    cseq  = seq[0:l]
    cname = "{}-len_{:02d}".format(_bname,l)
    pname = "{}-len_{:02d}".format(_bname,l-1)
    sfile = sfiles[e]

    if os.path.exists(pname+'.bar'):
      with open(pname+'.bar', 'r') as oldbar, \
        open(pname+'.map', 'w') as mapfile:
        for e, line in enumerate(oldbar):
          if e == 0: continue
          cols = line.strip().split()
          mapfile.write(cols[1] + '.' + "\n")
      mfile = pname+'.map'
    else :
      mfile = ''

    # Make sure the first round for mapping is always recomputed
    # force = True if e == 1 else args.force

    [sfile, bfile, efile, rfile, psfile] = ril.sys_barriers(cname, cseq, sfile, 
        barriers=args.barriers,
        minh=args.b_minh,
        maxn=args.b_maxn,
        k0=args.k0,
        temp=args.temperature,
        noLP=args.noLP,
        moves='single-base-pair',
        gzip=True,
        rates=True,
        bsize=False,
        saddle=False,
        mfile=mfile,
        force=args.force,
        verb=args.verbose)
    bfiles.append([bfile, efile, rfile, psfile])
  return bfiles

def barmap_subopts(_sname, seq, args):
  """ Print all suboptimal structures, return list of files """
  sfiles = []
  for l in range(args.start, args.stop+1):
    cseq  = seq[0:l]
    cname = "{}-len_{}".format(_sname,l)
  
    csfile = ril.sys_suboptimals(cname, cseq, 
        RNAsubopt=args.RNAsubopt,
        ener=args.s_ener, 
        temp=args.temperature,
        noLP=args.noLP,
        opts=args.spatch,
        verb=args.verbose, 
        force=args.force)
    sfiles.append(csfile)
  return sfiles

def set_p0(bfile, l, lastlines, curlmin, newlmin, cutoff, verb):
  '''
  remap densities from :lastlines: using the information in plist
  '''

  lminmap = c.defaultdict(int)
  for x,y in zip(curlmin, newlmin): 
    lminmap[x]=y

  p0dict = c.defaultdict(float)
  for (i,pop) in enumerate(lastlines[0].split()):
    if i == 0 :
      time = float(pop)
    else :
      p0dict[lminmap[i]]+=float(pop)
      if float(pop) > cutoff:
        print "{:3d} {:3d} {:f} {:s} {:d} => {:d}".format(
            l, i, float(pop), get_structure(bfile, i), i, lminmap[i])
        if lminmap[i] == 0:
          sys.exit('lost significant population!')

  p0 = [] 
  p0sum = 0.0
  for (x,y) in p0dict.items():
    if x == 0 : continue
    if y > cutoff :
      p0.append("{:d}={:f}".format(x, y))
      p0sum += y

  if verb : 
    print >> sys.stderr, "Total population {:.3f}\n".format(p0sum)
  return p0

def get_structure(bfile, idx):
  ss=''
  with open(bfile, 'r') as bar:
    for n, line in enumerate(bar):
      if n == idx : 
        ss = line.strip().split()[1]
        break
  return ss

def get_mapping_dict(oldbar, minfo):
  '''
  * parse the old barfile
      o = [[idx, struct],...]
  * parse the current barfile
      c = [[idx, struct],...]
  * dictionary[old_idx]=[new_idx,dist]
  '''
  mapper = c.defaultdict(int) 
  with open(oldbar, 'r') as old, \
      open(minfo, 'r') as info:

    mapinfo = []
    for n, line in enumerate(info):
      if re.match('[\.\(\)]+', line.strip()):
        gstr, sptidx, energy, fmin, fminT, gmin, gminT = line.strip().split()
        mapinfo.append([gminT, gstr])
      elif re.match('not in hash', line.strip()):
        print >> sys.stderr, "structure not in hash"
        mapinfo.append([0, ''])
      elif re.match('not yet assigned', line.strip()):
        print >> sys.stderr, "structure not yet assigned"
        mapinfo.append([0, ''])

    for n, line in enumerate(old):
      if n == 0 : continue
      else : # old idx = enumerate n
        [old_idx, old_min] = line.strip().split()[0:2]
        if n > len(mapinfo):
          sys.exit('ERROR: To many lines in barfile: ' \
              + minfo + ' missing mapping information!')
        mapper[int(old_idx)]=int(mapinfo[n-1][0])

  return mapper

def pathlist(mapdata):
  '''
  Translate a list of dictionaries with mapping information into
  a list of trajectories for lmin-transitions

  :return: plist
  '''

  # initialize path-lists with first lmins
  plist = [[i+1] for i in range(len(mapdata[0]))]
  
  # append the mapping for each current minimum
  for m, mdict in enumerate(mapdata):
    seen = set()
    for path in plist:
      last = path[-1]
      #print path, "last: ", last, mdict[last]
      path.append(mdict[last])
      seen.add(mdict[last])

    #for path in plist: print "update: ", path

    # add new pathlists for lmins that appear in 
    # the next landscape for the first time
    if m == len(mapdata)-1 : break
    for newmin in sorted(mapdata[m+1], key=mapdata[m+1].get, reverse=False):
      if newmin not in seen:
        newpath = [0 for i in range(len(plist[0])-1)]
        newpath.append(newmin)
        #print "Adding new path: ", newpath
        plist.append((newpath))

  return sorted(plist, key=lambda m:m[-1])

def get_barmap_args():
  """ A collection of arguments that are used by BarMap """
  parser = argparse.ArgumentParser()

  ril.argparse_add_arguments(parser, 
      RNAsubopt=True,
      barriers=True,
      treekin=True,
      noLP=True, temperature=True, 
      tmpdir=True, name=True, force=True, verbose=True,
      start=True, stop=True, k0=True, tX=True, cutoff=True)

  parser.add_argument("--plot_cutoff", type=float, default=0.02)
  parser.add_argument("--plot_title", default='')
  parser.add_argument("--noplot", action="store_true", 
      help="Do not plot results") 
  parser.add_argument("--plot_linlog", action="store_true",
      help="Divide x-axis into lin and log at transcription stop")

  ##############
  # Depricated:
  #TODO: need to think whether this is the best way to do it... hidden, for now!
  parser.add_argument("--spatch", 
      default='pylands_spatch.py', action = 'store',
      help=argparse.SUPPRESS)
      #help="Specify a script to postprocess suboptimals")
  parser.add_argument("--theo", action="store_true",
      help="Use *hidden* spatch script to include theophylline contributions")
  #parser.add_argument("--s_stmp", default='/tmp', action = 'store',
  #    help="Specify path to a temporary sort directory for unix sort")

  return parser.parse_args()

def main():
  """ BarMap-v2.0 -- cotransriptional folding
    Dependencies: RNAsubopt, barriers, treekin
    
    The implementation is split into 4 steps ...
  """

  ### parse input & adjust arguments ###
  args = get_barmap_args()

  # Read Input & Update Arguments
  name, seq = ril.parse_vienna_stdin(sys.stdin)

  # One name, just to be clear ...
  (args.name, name) = (args.name, args.name) if args.name else (name, name)

  if args.stop == 0 : 
    args.stop = len(seq)
  else :
    seq = seq[:args.stop]

  if args.verbose: 
    print "# Input: {:s} {:s}".format(name, seq)
    
  if args.s_ener is None :
    args.s_ener, args.s_maxn = ril.sys_subopt_range(seq, 
        nos=args.s_maxn, maxe=args.s_maxe, verb=(args.verbose>0))
    if args.verbose:
      print "# Energyrange {:.2f} computes {:d} sequences".format(
          args.s_ener, args.s_maxn)
  elif args.verbose : 
    args.s_ener, args.s_maxn = ril.sys_subopt_range(seq, 
        nos=0, maxe=args.s_ener, verb=False)
    print "# Energyrange {:.2f} computes {:d} sequences".format(
        args.s_ener, args.s_maxn)

  # Spatch-Hack to include theophylline binding
  if which(args.spatch) is None :
    print >> sys.stderr, args.spatch, "is not executable"
    print """ 
    You need to make sure that the executables of pylands are in your
    $PATH environment.
    """
    raise RuntimeError('Could not find executable')

  spatch = ['|', args.spatch]
  args.spatch = ''
  if args.theo:
    spatch.append('--theo')
    #spatch.append('-v')
    name += '_theo'; args.name = name
    args.spatch=spatch

  if not args.tmpdir :
    args.tmpdir = 'BarMap_'+args.name

  if not os.path.exists(args.tmpdir):
    os.makedirs(args.tmpdir)

  print """# Starting with BarMap computations ... """
  #TODO: assign names here in the main loop

  print """# writing RNAsubopt files ... """
  sname = "{}/{}-ener_{:.2f}".format(args.tmpdir, args.name, args.s_ener)
  #if args.circ: myfile += '_circ'
  if args.noLP: sname += '_noLP'
  sfiles = barmap_subopts(sname, seq, args)

  print """# writing barriers files ... """
  bname = "{}-minh_{}-maxn_{}-k0_{}".format(sname, 
      args.b_minh, args.b_maxn, args.k0)
  bfiles = barmap_barriers(bname, seq, sfiles, args)

  print """# writing/parsing mapping information ... """
  plist = barmap_mapping(bname, seq, args)

  print """# simulations using treekin ... """
  tfiles = barmap_treekin(bname, seq, bfiles, plist, args)

  if not args.noplot :
    print """# merging treekin results into BarMap plot ... """
    plotfile = barmap_plot(name, seq, tfiles, plist, args)
    print "# Your results have been plotted in the file: {}".format(plotfile) 

  print "# Thank you for using BarMap b(^.^)d"
  return

if __name__ == '__main__':
  main()

