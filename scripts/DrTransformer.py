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

def add_transition_edges(CG, saddles, fullseq, s1, s2, ts=None, 
    maxdG=10.00, _k0 = 2e5, _fpath = 10, _RT = 0.61632077549999997):
  """ compute the transition rates between two structures: s1 <-> s2, 
    where s2 is always the new, energetically better structure.
  """
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

  valid = (ts is not None or saddleE-CG.node[s1]['energy'] <= maxdG)

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

def dump_conformation_graph(CG, seq, name, verb=False) :
  """ Make a barriers + rates output from the current conformation graph """
  bfile = name+'.bar'
  rfile = name+'.rts'
  p0 = []

  sorted_nodes = sorted(CG.nodes(data=True), 
      key=lambda x: x[1]['energy'], reverse=False)
  with open(bfile, 'w') as bar, open(rfile, 'w') as rts :
    bar.write("     {}\n".format(seq))
    for e, (ni, data) in enumerate(sorted_nodes) :
      bar.write("{:4d} {} {:6.2f}\n".format(e+1, ni[:len(seq)], data['energy']))
      if verb :
        print "{:4d} {} {:6.2f} {:6.4f}".format(
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

def get_stats_and_update_occupancy(CG, sorted_nodes, tfile) :
  """
    Update the occupancy in the Graph and the total simulation time
  """

  # http://www.regular-expressions.info/floatingpoint.html
  reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')

  lastlines = s.check_output(['tail', '-2', tfile]).strip().split("\n")
  if not reg_flt.match(lastlines[0]):
    sys.exit('over and out')
  else :
    time = float(lastlines[0].split()[0])
    iterations = int(lastlines[-1].split()[-1])
    for e, occu in enumerate(lastlines[0].split()[1:]) :
      ss = sorted_nodes[e][0]
      CG.node[ss]['occupancy'] = float(occu)
  return time, iterations

def graph_pruning(CG, sorted_nodes, cutoff, saddles, fullseq) :
  """ Delete nodes or report them as still reachable """
  deleted_nodes = 0
  still_reachables = 0

  for ni, data in reversed(sorted_nodes) :
    #print ni, data
    if data['occupancy'] < cutoff :
      nbrs = sorted(CG.successors(ni), 
          key=lambda x: CG.node[x]['energy'], reverse=False)
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

def expand_graph(CG, seq, saddles, 
    maxdG = None, 
    _cutoff = 0.01, 
    mfe_only=False,
    verb=False):
  """ Expand the graph ...

  :return: Number of new nodes

  """
  csid = CG.graph['seqid']
  fseq = CG.graph['full_sequence']
  tlen = CG.graph['transcript_length']

  # Add MFE
  ss, mfe = RNA.fold(seq)
  future = '.' * (len(fseq)-tlen)
  ss = ss + future
  #print >> sys.stderr, "{}\n{} {:6.2f}".format(seq, ss, mfe)

  if CG.has_node(ss) :
    en = CG.node[ss]['energy']
  else :
    en = round(RNA.energy_of_structure(fseq, ss, 0), 2)
    if nx.number_of_nodes(CG) == 0 :
      CG.add_node(ss, energy=en, occupancy=1.0, identity=CG.graph['seqid'])
      CG.graph['seqid'] += 1
    else :
      for ni in CG.nodes() :
        if CG.has_node(ss):
          #print """ node exists """
          add_transition_edges(CG, saddles, fseq, ni, ss, maxdG=maxdG)
        elif add_transition_edges(CG, saddles, fseq, ni, ss, maxdG=maxdG) :
          CG.node[ss]['energy'] = en
          CG.node[ss]['occupancy'] = 0.0
          CG.node[ss]['identity'] = CG.graph['seqid']
          CG.graph['seqid'] += 1

  if not CG.has_node(ss) :
    print "# WARNING: ", ss, "[mfe secondary structure could not be connected]"

  if mfe_only is True :
    pass
  else :
    for ni, data in CG.nodes_iter(data=True):
      en  = data['energy']
      occ = data['occupancy']
      if occ < _cutoff : continue

      ss = ni[0:len(seq)]

      opened = open_breathing_helices(seq, ss)
      for onbr in opened :
        nbr = fold_exterior_loop(seq, onbr)
        future = '.' * (len(ni) - len(nbr))
        nbr += future

        if ni == nbr or CG.has_edge(ni, nbr):
          continue

        if CG.has_node(nbr):
          """ node exists """
          add_transition_edges(CG, saddles, fseq, ni, nbr, maxdG=maxdG)
        elif add_transition_edges(CG, saddles, fseq, ni, nbr, maxdG=maxdG) :
          """ node does not exist """
          enbr = round(RNA.energy_of_structure(fseq, nbr, 0), 2)
          CG.node[nbr]['energy'] = enbr
          CG.node[nbr]['occupancy'] = 0.0
          CG.node[nbr]['identity'] = CG.graph['seqid']
          CG.graph['seqid'] += 1
        else :
          print "# WARNING: Could not add transition edge!"
  return CG.graph['seqid']-csid

def open_breathing_helices(seq, ss, free=3):
  """ """
  nbrs = set()
  pt = nal.make_pair_table(ss, base=0)

  # mutable secondary structure 
  nbr = list(ss)

  rec_fill_nbrs(nbrs, ss, nbr, pt, (0, len(ss)), free)

  nbrs.add(''.join(nbr))

  return nbrs

def rec_fill_nbrs(nbrs, ss, mb, pt, (n, m), free):
  """ 
    TODO: Test function, but looks good
  """
  skip = 0 # fast forward in case we have deleted stuff
  for i in range(n, m) :
    j = pt[i]
    if j == -1 : continue
    if i < skip: continue

    nb = list(ss)
    [o,l] = [0,0]
    [p,q] = [i,j]

    add = True
    while p < q and (l == 0 or o < free):
      if pt[p] != q or p != pt[q] :
        """ this is a multiloop """
        # i,j = 1, len(pt)
        rec_fill_nbrs(nbrs, ''.join(nb), mb, pt, (p,q), free-o)
        add = False
        break
      else :
        # remove the base-pairs
        pt[p] = pt[q] = -1
        nb[p] = nb[q] = '.'
        mb[p] = mb[q] = '.'
        o += 2 # one base-pair deleted, two bases freed

      l = 0 # reset interior-loop size
      while (p < q and pt[p+1] == -1):
        p += 1
        l += 1
      p += 1
      while (p < q and pt[q-1] == -1):
        q -= 1
        l += 1
      q -= 1
      o += l

    if add :
      nbrs.add(''.join(nb))
    skip = j+1

  return 

def fold_exterior_loop(seq, con, fast=1):
  """ Constrained folding, if fast, only fold the unconstrained exterior loop
  of a structure, leave the rest as is. """

  if fast :
    spacer = 'NNN'
    pt = nal.make_pair_table(con, base=0)
    ext = ''

    # shrink the sequcnes
    skip = 0
    for i, j in enumerate(pt):
      if i < skip : continue
      if j == -1 : 
        ext += seq[i]
      else :
        ext += spacer
        skip = j+1
    css, cfe = RNA.fold(ext)
    
    # replace characters in constraint
    c, skip = 0, 0
    for i, j in enumerate(pt):
      if i < skip : continue
      if j == -1 : 
        con = con[:i] + css[c] + con[i+1:]
        c += 1
      else :
        c += len(spacer)
        skip = j+1
    ss = con

  else :
    # Force copy of string for ViennaRNA swig interface bug
    tmp = (con + '.')[:-1]
    RNA.cvar.fold_constrained = 1
    ss, mfe = RNA.fold(seq, tmp)
    RNA.cvar.fold_constrained = 0

  return ss

def talk_to_DrForna(CG, seq, sorted_nodes, tfile, _time, _cutoff, 
    repl=None, _outfile=None) :
  """  """
  # http://www.regular-expressions.info/floatingpoint.html
  reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')

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
            if repl and prevcourse and prevcourse[e+1] > 0 \
                and line is not tknlines[-2]:
              y1 = min(prevcourse[e+1], course[e+1])
              y2 = max(prevcourse[e+1], course[e+1])
              dy = math.log(y2) - math.log(y1)
              if dy < repl : continue
            ss = sorted_nodes[e][0]

            line = "{} {} {} {:s} {:6.2f}".format(CG.node[ss]['identity'], 
                _time + time, occu, ss[:len(seq)], CG.node[ss]['energy'])

            if _outfile :
              with open(_outfile, 'a') as outf:
                outf.write(line + '\n')
            else :
              print line
        prevcourse = course
  return 

def get_drtrafo_args():
  """ A collection of arguments that are used by DrTransformer """
  parser = argparse.ArgumentParser(
      #formatter_class=argparse.RawTextHelpFormatter,
      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
      #formatter_class=argparse.MetavarTypeHelpFormatter,
      description='echo sequence | %(prog)s [options]')

  parser.add_argument("--findpath", type = int, default = 10, metavar='<int>',
      help="Specify search width for *findpath* heuristic") 
  parser.add_argument("--minrate", type = float, default = 1e-10, 
      metavar='<flt>',
      help="Specify minmum rate to accept a new neighbor")
  parser.add_argument("--cutoff", type=float, default=0.01, metavar='<flt>',
      help="Cutoff for population transfer")

  parser.add_argument("--tmpdir", default='DrTrafo_Files', action = 'store', 
      metavar='<str>',
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
  parser.add_argument("-T","--temperature", type=float, default=37.0, 
      metavar='<flt>',
      help="Set the temperature in Celsius for ViennaRNA computations")

  parser.add_argument("--treekin", default='treekin', action = 'store',
      help="Specify path to your *treekin* executable")
  parser.add_argument("--repl", type=float, default=None,
      help="Logarithmically reduce output size")
  #parser.add_argument("--k0", type=float, default=2e5,
  #    help="Arrhenius rate prefactor")
  parser.add_argument("--t0", type=float, default=0.0,
      help="First time point of the printed time-course")
  parser.add_argument("--ti", type=float, default=1.02,
      help="Output-time increment of solver (t1 * ti = t2)")
  parser.add_argument("--t8", type=float, default=0.02,
      help="Transcription speed in seconds per nucleotide")
  parser.add_argument("--tX", type=float, default=60,
      help="Simulation time after transcription")

  #parser.add_argument("--noplot", action="store_true", 
  #    help="Do not plot results") 
  parser.add_argument("--plot_cutoff", type=float, default=0.02)
  #parser.add_argument("--plot_title", default='')
  #parser.add_argument("--plot_linlog", action="store_true",
  #    help="Divide x-axis into lin and log at transcription stop")

  return parser.parse_args()

def main():
  """ DrTransformer - cotranscriptional folding """
  args = get_drtrafo_args()

  (name, fullseq) = nal.parse_vienna_stdin()

  # Adjust arguments, prepare simulation
  if args.name == '' : 
    args.name = name
  else :
    name = args.name

  _outfile = name + '.drf'
  _RT= 0.61632077549999997
  if args.temperature != 37.0 :
    _RT = (_RT/37.0)*args.temperature
    print >> sys.stderr, 'WARNING: temperature option not fully supported'
  _last_only  = 0, # Print only the last simulation
  _output     = 'DrForna'#, 'BarMap'

  if args.stop == 0 : 
    args.stop = len(fullseq)+1
  else :
    fullseq = fullseq[0:args.stop-1]

  if not os.path.exists(args.tmpdir):
    os.makedirs(args.tmpdir)

  # Minrate specifies the lowest accepted rate for simulations (sec^-1)
  # it can be directly converted into a activation energy that results this rate
  maxdG = -_RT * math.log(args.minrate)
  #print args.minrate, '=>', maxdG


  # Start with DrTransformer
  RNA.cvar.noLonelyPairs = 1
  RNA.cvar.temperature = args.temperature

  if _output == 'DrForna':
    with open(_outfile, 'w') as outf :
      outf.write("id time conc struct energy\n")
  if args.verbose :
    print "T: i,j,pop,struct,energy"

  # initialize a directed conformation graph
  CG = nx.DiGraph(
      full_sequence=fullseq, 
      transcript_length=None,
      total_time=0,
      seqid=0)
  saddles = c.defaultdict()

  for tlen in range(args.start, args.stop) :
    CG.graph['transcript_length']=tlen
    seq = fullseq[0:tlen]

    # Expand Graph
    nnn = expand_graph(CG, seq, saddles, maxdG=maxdG, 
        _cutoff = args.cutoff, mfe_only=True, verb=args.verbose)

    # Report before simulation (?)
    # print nnn, "new nodes"

    # Simulate
    _fname = args.tmpdir+'/'+name+'-'+str(tlen)
    _t8 = args.tX if tlen == args.stop-1 else args.t8

    # produce input for treekin simulation
    [bfile, rfile, p0, nlist] = dump_conformation_graph(CG, seq, _fname, 
        verb=args.verbose)

    dn,sr = 0,0
    if len(nlist) == 1 :
      # Fake Results for DrForna
      CG.graph['total_time'] += _t8
      if _output == 'DrForna' :
        ss = nlist[0][0]
        line = "{} {} {} {:s} {:6.2f}".format(CG.node[ss]['identity'], 
            CG.graph['total_time'], 1.0, ss[:len(seq)], CG.node[ss]['energy'])
        with open(_outfile, 'a') as outf:
          outf.write(line + '\n')
    else :
      # - Simulate with treekin
      try:
        tfile = nal.sys_treekin(_fname, seq, bfile, rfile, 
            treekin=args.treekin, p0=p0, t0=args.t0, ti=args.ti, t8=_t8, 
            useplusI=True, force=True, verb=False)
      except RuntimeError:
        try :
          tfile = nal.sys_treekin(_fname, seq, bfile, rfile, 
              treekin=args.treekin, p0=p0, t0=args.t0, ti=args.ti, t8=_t8, 
              useplusI=False, force=True, verb=True)
        except RuntimeError:
          print "Abort after", tlen, "nucleotides:", \
              "treekin cannot find a solution, sorry"
          raise SystemExit

      # Get Results
      old_time = CG.graph['total_time']
      time_inc, iterations = get_stats_and_update_occupancy(CG, nlist, tfile)
      if _output == 'DrForna' :
        talk_to_DrForna(CG, seq, nlist, tfile, old_time,
          args.plot_cutoff, repl=args.repl, _outfile=_outfile)
      CG.graph['total_time'] += time_inc
      
      # Prune
      dn,sr = graph_pruning(CG, nlist, args.cutoff, saddles, fullseq)

    if args.verbose :
      print "# Deleted {} nodes, {} still reachable.".format(dn, sr)

  return

if __name__ == '__main__':
  main()

