#!/usr/bin/env python

"""
  Code written by Stefan Badelt (stef@tbi.univie.ac.at)

  vim-config = set: ts=2 et sw=2 sts=2
"""

import re
import os
import sys
import math
import argparse
import networkx as nx
import subprocess as s
import collections as c

import rnaworld as nal
import RNA

def add_transition_edges(CG, saddles, fullseq, s1, s2, ts=None) :
  """ compute the transition rates between two structures: s1 <-> s2, 
    where s2 is always the new, energetically better structure.
  """
  _RT    = 0.61632077549999997 # unit?
  _maxS  = 10.00  # kcal/mol
  _k0    = 2e5    # sec^-1
  _fpath = 10

  saddleE = None
  if (s1, s2) in saddles :
    saddleE = saddles[(s1,s2)]
  else :
    saddleE = float(RNA.find_saddle(fullseq, s1, s2, _fpath))/100

  # Minimum between direct and in-direct path barriers
  if ts : # then we know that the indirect path has to be part of saddles
    tsE1 = saddles[(s1,ts)] if (s1,ts) in saddles else saddles[(ts,s1)]
    tsE2 = saddles[(s2,ts)] if (s2,ts) in saddles else saddles[(ts,s2)]
    tsE = max(tsE1, tsE2)
    saddleE = min(tsE, saddleE)

  saddles[(s1,s2)] = saddleE

  valid = (ts is not None or saddleE <= _maxS)

  if valid :
    e1 = CG.node[s1]['energy']
    e2 = round(RNA.energy_of_structure(fullseq, s2, 0), 2)

    # Energy barrier
    dG_1s = round(saddleE-e1, 2)
    dG_2s = round(saddleE-e2, 2)

    # Metropolis Rule
    k_12 = _k0 * math.e**(-dG_1s/_RT) if dG_1s > 0 else _k0
    k_21 = _k0 * math.e**(-dG_2s/_RT) if dG_2s > 0 else _k0

    CG.add_weighted_edges_from([(s1, s2, k_12)])
    CG.add_weighted_edges_from([(s2, s1, k_21)])

  #print "#Added Edge:", s1, s2, "({}, {:g}, {:g})".format(valid, k_12, k_21)

  return valid

def dump_conformation_graph(CG, seq, name) :
  """ Make a barriers + rates output from the current conformation graph """
  bfile = name+'.bar'
  rfile = name+'.rts'
  p0 = []

  sorted_nodes = sorted(CG.nodes(data=True), key=lambda x: x[1]['energy'], reverse=False)
  with open(bfile, 'w') as bar, open(rfile, 'w') as rts :
    bar.write("     {}\n".format(seq))
    for e, (ni, data) in enumerate(sorted_nodes) :
      bar.write("{:4d} {} {:6.2f}\n".format(e+1, ni[:len(seq)], data['energy']))
      print >> sys.stderr, "{:4d} {} {:6.2f} {:6.4f}".format(
          e+1, ni[:len(seq)], data['energy'], data['occupancy'])
      if data['occupancy'] > 0 :
        p0.append("{}={}".format(e+1,data['occupancy']))
      rates = []
      for (nj, jdata) in sorted_nodes :
        if CG.has_edge(ni,nj) :
          rates.append(CG[ni][nj]['weight'])
        else :
          rates.append(0)
      line = "".join(map("{:10.4g}".format, rates))
      rts.write("{}\n".format(line))

  return [bfile, rfile, p0, sorted_nodes]

def update_occupancy(CG, sorted_nodes, tfile) :
  """
    Update the occupancy in the Graph and the total simulation time
  """

  # http://www.regular-expressions.info/floatingpoint.html
  reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')

  lastlines = s.check_output(['tail', '-2', tfile]).strip().split("\n")
  if not reg_flt.match(lastlines[0]):
    sys.exit('over and out')
  else :
    #time = float(lastlines[0].split()[0])
    #iterations = int(lastlines[-1].split()[-1])
    for e, occu in enumerate(lastlines[0].split()[1:]) :
      ss = sorted_nodes[e][0]
      CG.node[ss]['occupancy'] = float(occu)
  return

def graph_pruning(CG, sorted_nodes, cutoff, saddles, fullseq) :
  """ Delete nodes or report them as still reachable """
  deleted_nodes = 0
  still_reachables = 0

  for ni, data in reversed(sorted_nodes) :
    #print ni, data
    if data['occupancy'] < cutoff :
      nbrs = sorted(CG.successors(ni), key=lambda x: CG.node[x]['energy'], reverse=False)
      best, been = nbrs[0], CG.node[nbrs[0]]['energy']

      if been > data['energy'] :
        still_reachables += 1
        #print "still reachable", best, been
        continue

      for e, nbr in enumerate(nbrs[1:]) :
        check = add_transition_edges(CG, saddles, fullseq, nbr, best, ni)
        if not check :
          print e, nbr, best, ni, check
          sys.exit('over and out')

      CG.remove_node(ni)
      deleted_nodes += 1

  return deleted_nodes, still_reachables

def expand_graph(CG, seq, ss, _cutoff = 0.01, verb=False):
  """ Expand the graph ...

  """
  for ni, data in CG.nodes_iter(data=True):
    en  = data['energy']
    occ = data['occupancy']
    if occ < _cutoff : continue

    ss = ni[0:len(seq)]

    opened = open_breathing_helices(seq, ss)
    for onbr in opened :
      nbr = fold_exterior_loop(seq, onbr)
      future = '.' * (ni-nbr)
      nbr += future
      if ni == nbr or CG.has_edge(ni, nbr):
        continue

      if CG.has_node(nbr):
        add_transition_edges(CG, saddles, fullseq, ni, nbr)
      elif add_transition_edges(CG, saddles, fullseq, ni, nbr) :
        enbr = round(RNA.energy_of_structure(fullseq, nbr, 0), 2)
        CG.node[nbr]['energy'] = enbr
        CG.node[nbr]['occupancy'] = 0.0
        CG.node[nbr]['identity'] = _id
        _id += 1

  return 

def open_breathing_helices(seq, ss):

  return []

def fold_exterior_loop(seq, onbr):

  return onbr

def talk_to_DrForna(CG, seq, sorted_nodes, tfile, _time, _cutoff, repl=None) :
  """ TODO: this should not return time, it should just talk! """
  # http://www.regular-expressions.info/floatingpoint.html
  reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')

  if len(sorted_nodes) == 1 :
    time = 0
    ss = sorted_nodes[0][0]
    print "{} {} {} {:s} {:6.2f}".format(CG.node[ss]['identity'], 
        _time, 1.0, ss[:len(seq)], CG.node[ss]['energy'])
  else :
    prevcourse = []
    with open(tfile) as tkn :
      # this is nasty, but used to check if we're at the last line
      tknlines = tkn.readlines()
      for line in tknlines:
        if reg_flt.match(line) :
          course = map(float, line.strip().split())
          time = course[0]

          for e, occu in enumerate(course[1:]) :
            # is it above visibility threshold?
            if occu > _cutoff :
              # can we actually skip the point?
              if repl and prevcourse and prevcourse[e+1] > 0 and line is not tknlines[-2]:
                y1 = min(prevcourse[e+1], course[e+1])
                y2 = max(prevcourse[e+1], course[e+1])
                dy = math.log(y2) - math.log(y1)
                if dy < repl : continue
              ss = sorted_nodes[e][0]
              print "{} {} {} {:s} {:6.2f}".format(CG.node[ss]['identity'], 
                  _time + time, occu, ss[:len(seq)], CG.node[ss]['energy'])

          prevcourse = course

  return time + _time

def get_drtrafo_args():
  """ A collection of arguments that are used by DrTransformer """
  parser = argparse.ArgumentParser(
      #formatter_class=argparse.RawTextHelpFormatter,
      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
      #formatter_class=argparse.MetavarTypeHelpFormatter,
      description='echo sequence | %(prog)s [options]')

  parser.add_argument("--findpath", type = int, default = 10, metavar='<int>',
      help="Specify search width for *findpath* heuristic") 
  parser.add_argument("--minrate", type = float, default = 1e-10, metavar='<flt>',
      help="Specify minimum rate for accepting new neighbors")
  parser.add_argument("--cutoff", type=float, default=0.01, metavar='<flt>',
      help="Cutoff for population transfer")

  parser.add_argument("--tmpdir", default='DrTrafo_Files', action = 'store', metavar='<str>',
      help="Specify path for storing temporary output files")
  parser.add_argument("-v","--verbose", action="store_true",
      help="Track process by writing verbose output during calculations") 
  #parser.add_argument("-f","--force", action="store_true",
  #    help="Force to overwrite existing BarMap files") 
  parser.add_argument("-n","--name", default='', metavar='<string>',
      help="Name your output files, this option overwrites the fasta-header")
  parser.add_argument("--start", type=int, default=1, metavar='<int>',
      help="Start transcription at this nucleotide")
  parser.add_argument("--stop", type=int, default=0, metavar='<int>',
      help="Stop transcription at this nucleotide")
  parser.add_argument("-T","--temperature", type=float, default=37.0, metavar='<flt>',
      help="Set the temperature in Celsius for ViennaRNA computations")

  parser.add_argument("--treekin", default='treekin', action = 'store',
      help="Specify path to your *treekin* executable")
  parser.add_argument("--repl", type=float, default=None)
  parser.add_argument("--p0", nargs='+', default=['1=1'])
  parser.add_argument("--k0", type=float, default=2e5)
  parser.add_argument("--t0", type=float, default=0.0)
  parser.add_argument("--ti", type=float, default=1.02)
  parser.add_argument("--t8", type=float, default=0.02)
  parser.add_argument("--tX", type=float, default=60)

  parser.add_argument("--noplot", action="store_true", 
      help="Do not plot results") 
  parser.add_argument("--plot_cutoff", type=float, default=0.02)
  parser.add_argument("--plot_title", default='')
  parser.add_argument("--plot_linlog", action="store_true",
      help="Divide x-axis into lin and log at transcription stop")

  return parser.parse_args()

def main():
  """ DrTransformer - cotranscriptional folding """
  args = get_drtrafo_args()

  if not os.path.exists(args.tmpdir):
    os.makedirs(args.tmpdir)

  _RT= 0.61632077549999997
  _last_only  = 0, # Print only the last simulation
  _output     = 'DrForna', 'BarMap'
  _replot     = 0.01, # logarithmically reduce plot size
  
  (name, fullseq) = nal.parse_vienna_stdin()

  # initialize a directed conformation graph
  CG = nx.DiGraph()

  saddles = c.defaultdict()

  # init_simulation_map

  print "id time conc struct energy"
  print >> sys.stderr, "T: i,j,pop,struct,energy"

  RNA.cvar.noLonelyPairs = 1
  RNA.cvar.temperature = args.temperature

  if args.stop == 0 : 
    args.stop = len(fullseq)+1
  else :
    fullseq = fullseq[0:args.stop-1]

  _id = 0
  _total_time = 0

  for l in range(args.start, args.stop) :
    seq = fullseq[0:l]

    # Add MFE
    ss, mfe = RNA.fold(seq)
    future = '.' * (args.stop-l-1)
    ss = ss + future
    #print >> sys.stderr, "{}\n{} {:6.2f}".format(seq, ss, mfe)

    if CG.has_node(ss) :
      en = CG.node[ss]['energy']
    else :
      en = round(RNA.energy_of_structure(fullseq, ss, 0), 2)
      if nx.number_of_nodes(CG) == 0 :
        CG.add_node(ss, energy=en, occupancy=1.0, identity=_id)
        _id += 1
      else :
        for ni in nx.nodes(CG) :
          if add_transition_edges(CG, saddles, fullseq, ni, ss) :
            CG.node[ss]['energy'] = en
            CG.node[ss]['occupancy'] = 0.0
            CG.node[ss]['identity'] = _id
            _id += 1

    if not CG.has_node(ss) :
      print >> sys.stderr, ss, "[secondary structure could not be connected]"

    # Expand Graph
    # expand_graph(CG, seq, ss, verb=args.verbose)

    # Report before simulation (?)

    # Simulate
    _fname = args.tmpdir+'/'+name+'-'+str(l)
    _t8 = args.tX if l == args.stop-1 else args.t8

    [bfile, rfile, p0, nlist] = dump_conformation_graph(CG, seq, _fname)

    if len(nlist) == 1 :
      # Fake Results
      _total_time += _t8
      talk_to_DrForna(CG, seq, nlist, None, _total_time, args.plot_cutoff)
      print >> sys.stderr, "# Deleted {} nodes, {} still reachable.".format(0, 0)
      continue

    # - Simulate with treekin
    try:
      tfile = nal.sys_treekin(_fname, seq, bfile, rfile, 
          treekin=args.treekin, p0=p0, t0=args.t0, ti=args.ti, t8=_t8, 
          useplusI=True, force=True, verb=False)
    except RuntimeError:
      tfile = nal.sys_treekin(_fname, seq, bfile, rfile, 
          treekin=args.treekin, p0=p0, t0=args.t0, ti=args.ti, t8=_t8, 
          useplusI=False, force=True, verb=True)

    # Get Results
    update_occupancy(CG, nlist, tfile)
    _total_time = talk_to_DrForna(CG, seq, nlist, tfile, _total_time,
        args.plot_cutoff, repl=args.repl)
    
    # Prune
    dn,sr = graph_pruning(CG, nlist, args.cutoff, saddles, fullseq)
    print >> sys.stderr, "# Deleted {} nodes, {} still reachable.".format(dn, sr)

  return

if __name__ == '__main__':
  main()

