#!/usr/bin/env python

# written by Stefan Badelt (stef@tbi.univie.ac.at)
# vim-config = set: ts=2 et sw=2 sts=2

import re
import os
import sys
import math
import argparse
import networkx as nx
import subprocess as s
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import shutil
import tempfile
import contextlib
from struct import pack

import ribolands as ril
from ribolands.syswraps import SubprocessError, ExecError, check_version
from ribolands.crnwrapper import DiGraphSimulator
from ribolands.trafo import fold_exterior_loop

import RNA

@contextlib.contextmanager
def smart_open(filename=None, mode='w'):
  if filename and filename != '-':

    fh = open(filename, mode)
  else:
    fh = sys.stdout

  try:
    yield fh
  finally:
    if fh is not sys.stdout:
      fh.close()

def add_transition_edges(CG, saddles, args, s1, s2, ts=None): 
  """ compute the transition rates between two structures: s1 <-> s2, 
    where s2 is always the new, energetically better structure.
    NOTE: actually, that is not always the case if you use exterior_folding,
    then, s1 is always the conformation that is already present in the confomration
    graph...
  """
  maxdG=args.maxdG
  k0 = args.k0
  fpath = args.findpath_width
  _RT = args._RT
  fullseq = CG.graph['full_sequence']
  md = CG.graph['model_details']
  fc = CG.graph['fold_compound']

  saddleE = None
  if (s1, s2) in saddles :
    saddleE = saddles[(s1,s2)]
  else :
    saddleE = float(fc.path_findpath_saddle(s1, s2, fpath))/100
    #saddleE = float(RNA.find_saddle(fullseq, s1, s2, fpath))/100

  # Minimum between direct and in-direct path barriers
  if ts : # then we know that the indirect path has to be part of saddles
    tsE1 = saddles[(s1,ts)] #if (s1,ts) in saddles else saddles[(ts,s1)]
    tsE2 = saddles[(s2,ts)] #if (s2,ts) in saddles else saddles[(ts,s2)]
    tsE = max(tsE1, tsE2)
    saddleE = min(tsE, saddleE)

  saddles[(s1,s2)] = saddleE
  saddles[(s2,s1)] = saddleE

  valid = (ts is not None or saddleE-CG.node[s1]['energy'] <= maxdG)

  if valid :
    e1 = CG.node[s1]['energy']
    e2 = round(fc.eval_structure(s2), 2)
    #e2 = round(RNA.energy_of_structure(fullseq, s2, 0), 2)

    saddleE = max(saddleE, max(e1,e2)) # ensure saddle is not lower than s1, s2

    # Energy barrier
    dG_1s = saddleE-e1
    dG_2s = saddleE-e2

    # Metropolis Rule
    k_12 = k0 * math.exp(-dG_1s/_RT)
    k_21 = k0 * math.exp(-dG_2s/_RT)

    CG.add_weighted_edges_from([(s1, s2, k_12)])
    CG.add_weighted_edges_from([(s2, s1, k_21)])

    #print "#Added Edge:", s1, s2, "({}, {:g}, {:g})".format(valid, k_12, k_21)

  return valid

def dump_conformation_graph(CG, seq, saddles, name, logf=sys.stdout,verb=False) :
  """ Make a barriers + rates output from the current conformation graph.

  The printed files are input files for treekin. *.bar files contain the energy
  barriers to transition between local minima. Although every path eventually
  leads to the MFE structure, it can proceed via an energetically worse
  structure first. This is in contrast to files produced by `barriers`, where
  local minima are always *directly* connected to energetically better local
  minima.   
  """

  sorted_nodes = filter(lambda (n,d): d['active'], 
      sorted(CG.nodes(data=True), 
      key=lambda x: x[1]['energy'], reverse=False))

  barfile_nodes = filter(lambda d: CG.node[d]['active'], 
      sorted(CG.nodes(data=False), key=lambda x: CG.node[x]['energy'], reverse=False))

  if name :
    bfile = name+'.bar'
    rfile = name+'.rts'
    brfile = rfile+'.bin'
    p0 = []
    
    with open(bfile, 'w') as bar, open(rfile, 'w') as rts, open(brfile, 'w') as brts :
      bar.write("     {}\n".format(seq))
      brts.write(pack("i", len(sorted_nodes)))
      for e, (ni, data) in enumerate(sorted_nodes, 1) :

        # Calculate barrier heights to all other basins.
        nextmin = 0
        barrier = 0
        saddleE = None
        nMsE = set()
        for ee, be in enumerate(barfile_nodes, 1):
          if e == ee :
            continue
          if (ni, be) in saddles :
            sE = saddles[(ni,be)]
            assert sE == saddles[(be,ni)]
            nMsE.add((ee, sE))
          else :
            sE = None

        mystr = ' '.join(map(lambda(x,y):'({:3d} {:6.2f})'.format(x,y-data['energy']), 
            sorted(list(nMsE), key=lambda x:x[0])))

        bar.write("{:4d} {} {:6.2f} {}\n".format(e, ni[:len(seq)], data['energy'], 
          mystr))

        if verb :
          line = "{:4d} {:4d} {} {:6.2f} {:6.4f} (ID = {:d})\n".format(
              CG.graph['transcript_length'], e, ni[:len(seq)], 
              data['energy'], data['occupancy'], data['identity'])
          logf.write(line)

        if data['occupancy'] > 0 :
          p0.append("{}={}".format(e,data['occupancy']))
        trates = []
        rates = []
        for (nj, jdata) in sorted_nodes :
          if CG.has_edge(ni,nj) :
            rates.append(CG[ni][nj]['weight'])
            trates.append(CG[nj][ni]['weight'])
          else :
            rates.append(0)
            trates.append(0)
        line = "".join(map("{:10.4g}".format, rates))
        rts.write("{}\n".format(line))
        for r in trates:
          brts.write(pack("d", r))
  else :
    line = "Distribution of structures at the end:\n"
    for e, (ni, data) in enumerate(sorted_nodes, 1) :
      line += "LAST {:4d} {} {:6.2f} {:6.4f} (ID = {:d})\n".format(e, 
          ni[:len(seq)], data['energy'], data['occupancy'], data['identity'])
    logf.write(line)

    return 
  return [bfile, brfile, p0, sorted_nodes]

def get_stats_and_update_occupancy(CG, sorted_nodes, tfile) :
  """
    Update the occupancy in the Graph and the total simulation time
  """
  # http://www.regular-expressions.info/floatingpoint.html
  reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')

  lastlines = s.check_output(['tail', '-2', tfile]).strip().split("\n")
  if not reg_flt.match(lastlines[0]):
    raise ValueError('Cannot parse simulation output', tfile)
  else :
    time = float(lastlines[0].split()[0])
    iterations = int(lastlines[-1].split()[-1])
    tot_occ =sum(map(float, lastlines[0].split()[1:]))
    for e, occu in enumerate(lastlines[0].split()[1:]) :
      ss = sorted_nodes[e][0]
      CG.node[ss]['occupancy'] = float(occu)/tot_occ

  return time, iterations

def graph_pruning(CG, sorted_nodes, saddles, args) :
  """ Delete nodes or report them as still reachable. """
  cutoff = args.occupancy_cutoff

  deleted_nodes = 0
  still_reachables = 0

  for ni, data in reversed(sorted_nodes) :
    #print ni, data
    if data['occupancy'] < cutoff :

      nbrs = filter(lambda x: CG.node[x]['active'], sorted(CG.successors(ni), 
          key=lambda x: CG.node[x]['energy'], reverse=False))
      best, been = nbrs[0], CG.node[nbrs[0]]['energy']

      if been > data['energy'] :
        still_reachables += 1
        continue
      
      multibest = {best : been}
      #epsilon = 10 # kcal/mol
    
      (transfer, minbar) = (best, None)
      for e, nbr in enumerate(nbrs[1:]) :
        for mb in multibest.keys():
          always_true = add_transition_edges(CG, saddles, args, nbr, mb, ni)
          msE = saddles[(nbr, mb)]
          if minbar is None or minbar < (msE - CG.node[ni]['energy']):
            (transfer, minbar) = (nbr, msE - CG.node[ni]['energy'])
        #if CG.node[nbr]['energy'] <= been + epsilon:
        multibest[nbr] = CG.node[nbr]['energy']
        if always_true is False :
          raise ValueError('Did not add the transition edge!')

      if True: # Set to 'False' to keep all nodes
        CG.node[ni]['active']=False
        CG.node[transfer]['occupancy'] += CG.node[ni]['occupancy']
        CG.node[ni]['occupancy']=0.0
        deleted_nodes += 1

  return deleted_nodes, still_reachables

def expand_graph(CG, saddles, args, mode='default'):
  """ Find new neighbors and add them to the Conformation Graph

  The function is devided into two parts. 1) The current mfe structure
  is connected to all present structures, 2) The conformation graph is
  expanded using helix-breathing.

  :param CG: Conformation Graph (NetworkX)
  :param saddles: dictionary of all previous findpath runs
  :param args: commandline arguments and other global variables
    (using: cutoff, verbose)
  :param mode: choose from (1) mfe-only: only use current mfe as potential new
    neighbor (2) breathing-only: only use breathing neighborhood, (3) default:
    do both mfe and breathing.

  :return: Number of new nodes

  """
  cutoff= args.occupancy_cutoff
  verb  = args.verbose
  mfree = args.min_breathing

  csid = CG.graph['seqid']
  fseq = CG.graph['full_sequence']
  tlen = CG.graph['transcript_length']
  md = CG.graph['model_details']
  fc_full = CG.graph['fold_compound']
  seq = fseq[0:tlen]

  if mode not in ['default', 'mfe-only', 'breathing-only']:
    raise ValueError('unknown expansion mode')

  # Add MFE
  fc_tmp = RNA.fold_compound(seq, md)
  ss, mfe = fc_tmp.mfe()
  #ss, mfe = RNA.fold(seq)
  future = '.' * (len(fseq)-tlen)
  ss = ss + future
  #print >> sys.stderr, "{}\n{} {:6.2f}".format(seq, ss, mfe)

  regular_mode = True # NOTE: HACK! this is only here to produce any possible graph

  # If there is no node bec we are in the beginning, add the node,
  # otherwise, go through all nodes and try to add transition edges

  if nx.number_of_nodes(CG) == 0 :
    en = round(fc_full.eval_structure(ss), 2)
    #en = round(RNA.energy_of_structure(fseq, ss, 0), 2)
    CG.add_node(ss, energy=en, occupancy=1.0, 
        identity=CG.graph['seqid'], active=True)
    CG.graph['seqid'] += 1
  elif mode == 'default' or mode == 'mfe-only':
    # Try to connect MFE to every existing state
    for ni in CG.nodes() :
      if CG.node[ni]['active'] == False : continue
      if ni == ss or CG.has_edge(ni,ss) : continue

      if CG.has_node(ss): # from a previous iteration
        if add_transition_edges(CG, saddles, args, ni, ss):
          CG.node[ss]['active'] = True # in case it was there but inactive
      elif add_transition_edges(CG, saddles, args, ni, ss) :
        en = round(fc_full.eval_structure(ss), 2)
        #en = round(RNA.energy_of_structure(fseq, ss, 0), 2)
        CG.node[ss]['active'] = True
        CG.node[ss]['energy'] = en
        CG.node[ss]['occupancy'] = 0.0
        CG.node[ss]['identity'] = CG.graph['seqid']
        CG.graph['seqid'] += 1

  if mode == 'default' or mode == 'breathing-only':
    # Do the helix breathing graph expansion
    ext_moves = dict()
    # ext_moves[sequence] = [set((con,paren),...), structure]
    # sequence = exterior-loop sequence with NNN replacing constrained elements
    for ni, data in CG.nodes_iter(data=True):
      if data['active'] == False : continue
      en  = data['energy']
      occ = data['occupancy']
      if regular_mode and occ < cutoff : continue

      sss = ni[0:len(seq)]

      opened = open_breathing_helices(seq, sss, free=mfree)
      #print 'opened', opened
      for onbr in opened :
        nbr, ext = fold_exterior_loop(md, seq, onbr, ext_moves)
        future = '.' * (len(ni) - len(nbr))
        nbr += future

        if ni == nbr or CG.has_edge(ni, nbr):
          continue

        if CG.has_node(nbr):
          if add_transition_edges(CG, saddles, args, ni, nbr):
            CG.node[nbr]['active'] = True # in case it was there but inactive
        elif add_transition_edges(CG, saddles, args, ni, nbr) :
          enbr = round(fc_full.eval_structure(nbr), 2)
          #enbr = round(RNA.energy_of_structure(fseq, nbr, 0), 2)
          CG.node[nbr]['energy'] = enbr
          CG.node[nbr]['active'] = True
          CG.node[nbr]['occupancy'] = 0.0
          CG.node[nbr]['identity'] = CG.graph['seqid']
          CG.graph['seqid'] += 1
        else :
          """# WARNING: Could not add transition edge!"""

        # TODO: figure out if this is exact before merging into master
        if ext_moves[ext][0] :
          for (parent, child) in ext_moves[ext][0] :
            assert parent != ni # Parents may never be the same
            # Children can be the same, if the change is within the helix-breathing-move-set
            # p1 .((((((((.(((.......))).))).)))))......................
            # p2 ..(((((((.(((.......))).))).)))).......................
            # c1 .((((((((.(((.......))).))).)))))......(((......)))....
            # c2 .((((((((.(((.......))).))).)))))......(((......)))....
            if child == nbr : continue
            if CG.has_edge(parent, ni) :
              if CG.has_node(child) and CG.has_node(nbr):
                if CG.has_edge(nbr, child): 
                  continue
                # Calculate saddleE from saddleE of parents?
                sP = saddles[(parent, ni)] if (parent, ni) in saddles else saddles[(ni, parent)]
                sC1 = round(CG.node[child]['energy'] + sP - CG.node[parent]['energy'], 2)
                sC2 = round(CG.node[nbr]['energy'] + sP - CG.node[ni]['energy'], 2)

                if sC1 == sC2 and (nbr, child) not in saddles:
                    saddles[(nbr,child)] = sC1
                    saddles[(child,nbr)] = sC1
                #else :
                #  print CG.node[child]['energy'], sP, CG.node[parent]['energy']
                #  print CG.node[nbr]['energy'], sP, CG.node[ni]['energy']

                if add_transition_edges(CG, saddles, args, nbr, child):
                  CG.node[nbr]['active'] = True # in case it was there but inactive
                  CG.node[child]['active'] = True # in case it was there but inactive
                #else :
                #  print ext
                #  print 'sad', saddles[(nbr,child)]
                #  print nbr, CG.node[nbr]['energy']
                #  print child, CG.node[child]['energy']
                #  print saddles[(nbr,child)], CG[nbr][child]['weight'], CG[nbr][child]['weight']
                #  raise Exception("didn't add an expected neighbor")
                  
        # Track the final structure, every new identical ext-change will be
        # connected, if the parents were connected.
        ext_moves[ext][0].add((ni, nbr))

  if not CG.has_node(ss) or CG.node[ss]['active'] is False:
    print "# WARNING: ", ss, "[mfe secondary structure not connected]"

  return CG.graph['seqid']-csid

def open_breathing_helices(seq, ss, free=6):
  """ open all breathable helices, i.e. those that share a base-pair
    with an exterior loop region 
  """
  nbrs = set()
  pt = ril.make_pair_table(ss, base=0)

  # mutable secondary structure 
  nbr = list(ss)

  rec_fill_nbrs(nbrs, ss, nbr, pt, (0, len(ss)), free)

  nbrs.add(''.join(nbr))

  return nbrs

def rec_fill_nbrs(nbrs, ss, mb, pt, (n, m), free):
  """ recursive helix opening
  TODO: Test function, but looks good

  :param nbrs: a set of all neighboring conformations
  :param ss: reference secondary structure
  :param mb: a mutable version of ss, which, after the final round will have
    all breathing helices opened
  :param pt: pair table (zero based)
  :param (n,m): the range of the pt under current investigation
  :param free: number of bases that should be freed

  :return: 
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

      # remove the base-pairs
      pt[p] = pt[q] = -1
      nb[p] = nb[q] = '.'
      mb[p] = mb[q] = '.'
      o += 2 # one base-pair deleted, two bases freed

      l = 0 # reset interior-loop size
      while (p < q and pt[p+1] == -1):
        [p, l] = [p+1, l+1]
      p += 1
      while (p < q and pt[q-1] == -1):
        [q, l] = [q-1, l+1]
      q -= 1
      o += l

    if add :
      nbrs.add(''.join(nb))
    skip = j+1

  return 

def talk_generator(CG, sorted_nodes, tfile):
  """ Generator function that yields the time course information from the
  treekin output file.
  
  :param CG: Conformation Graph (networkx)
  :param sorted_nodes: a list of nodes sorted by their energy
  :param tfile: current treekin-output file

  """
  # http://www.regular-expressions.info/floatingpoint.html
  reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')

  ttime = CG.graph['total_time']

  with open(tfile) as tkn:
    # this is nasty, but used to check if we're at the last line
    prevcourse = []
    tknlines = tkn.readlines()
    for line in tknlines:
      if reg_flt.match(line) :
        course = map(float, line.strip().split())
        time = course[0]

        for e, occu in enumerate(course[1:]) :
          # is it above visibility threshold?
          ss = sorted_nodes[e][0]

          yield CG.node[ss]['identity'], ttime+time, occu, \
              ss[:CG.graph['transcript_length']], CG.node[ss]['energy']
        prevcourse = course
  return 

def plot_xmgrace(all_in, args):
  head = """
@with line
@line on
@line loctype world
@line g0
@line linewidth .2
@line linestyle 1
@line color 7
@line arrow 0
@line arrow type 0
@line arrow length 1.000000
@line arrow layout 1.000000, 1.000000
@line def
"""
  with open(args.name + '.gr', 'w') as gfh :
    gfh.write(head)
    for e, course in enumerate(all_in) :
      t, o = zip(*course)
      for i in range(len(t)) :
        gfh.write("{:f} {:f}\n".format(t[i], o[i]))
      gfh.write("&\n")

  return 

def plot_simulation(all_in, args):
  stop = args.stop
  start= args.start

  t8 = args.t8
  lin_time = (stop-start)*float(t8)
  tX = lin_time + args.tX if (lin_time + args.tX) >= lin_time * 10 else lin_time * 10

  name  = args.name
  title = args.name

  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)

  ax.set_ylim([0,1.01])
  ax.set_xscale('linear')

  # Make the second part of the plot logarithmic
  ax.set_xlim((0, lin_time))
  divider = make_axes_locatable(ax)
  axLog = divider.append_axes("right", size=2.5, pad=0, sharey=ax)
  axLog.set_xscale('log')
  axLog.set_xlim((lin_time+0.00001, tX))
  axLog.set_ylim([0,1.01])
  axLog.yaxis.set_visible(False)


  for e, course in enumerate(all_in) :
    if course == [] : continue
    t, o = zip(*course)
    # Determine which lines are part of the legend:
    # like this, it is only those that are populated 
    # at the end of transcription and if they reach 
    # an occupancy of 10% or higher
    if t[-1] > lin_time:
      p, = ax.plot(t, o, '-', lw=1.5)
      L, = axLog.plot(t, o, '-', lw=1.5)
      if max(o) >= 0.1 :
        L.set_label("ID {:d}".format(e))
    else :
      p, = ax.plot(t, o, '-', lw=0.5)
      L, = axLog.plot(t, o, '-', lw=0.5)

  fig.set_size_inches(7,3)
  fig.text(0.5,0.95, title, ha='center', va='center')

  #for tlen in range(args.stop-args.start) :
  #  ax.axvline(x=tlen*t8, linewidth=0.01, color='black', linestyle='--')

  # """ Add ticks for 1 minute, 1 hour, 1 day, 1 year """
  axLog.axvline(x=60, linewidth=1, color='black', linestyle='--')
  axLog.axvline(x=3600, linewidth=1, color='black', linestyle='--')
  axLog.axvline(x=86400, linewidth=1, color='black', linestyle='--')
  axLog.axvline(x=31536000, linewidth=1, color='black', linestyle='--')
  plt.legend()

  ax.set_ylabel('occupancy [mol/l]', fontsize=11)
  ax.set_xlabel('time [seconds]', ha='center', va='center', fontsize=11)
  ax.xaxis.set_label_coords(.9, -0.15)

  #plt.show()
  pfile = name+'.pdf'
  plt.savefig(pfile, bbox_inches='tight')

  return

def add_drtrafo_args(parser):
  """ A collection of arguments that are used by DrTransformer """
  # Treekin parameters
  parser.add_argument('--version', action='version', version='%(prog)s ' + ril.__version__)
  parser.add_argument("--treekin", default='treekin', action = 'store',
      metavar='<str>', help="Path to the *treekin* executable.")
  parser.add_argument("--ti", type=float, default=1.2, metavar='<flt>',
      help="""Output-time increment of treekin solver (t1 * ti = t2).""")

  # Common parameters
  parser.add_argument("--t8", type=float, default=0.02, metavar='<flt>',
      help="Transcription speed [seconds per nucleotide].")
  parser.add_argument("--tX", type=float, default=60, metavar='<flt>',
      help="Post-transcriptional simulation time [seconds].")
  parser.add_argument("-T","--temperature", type=float, default=37.0,
      metavar='<flt>', help="The temperature for ViennaRNA computations.")

  # Advanced algorithmic parameters
  parser.add_argument("--occupancy-cutoff", type=float, default=0.01, 
      metavar='<flt>',
      help="Occupancy cutoff for secondary structure graph expansion.")
  parser.add_argument("--findpath-width", type = int, default = 20, 
      metavar='<int>',
      help="Search width for *findpath* heuristic.") 
  parser.add_argument("--min-rate", type = float, default = 1e-10, 
      metavar='<flt>',
      help="""Minmum rate to accept a new structure as neighboring
      conformation.""")
  parser.add_argument('--structure-search-mode', default = 'default',
      choices=('default', 'mfe-only', 'breathing-only'),
      help="""Specify one of three modes: *default*: find new secondary
      structures using both the current MFE structure and breathing neighbors.
      *mfe-only*: only find the current MFE structure at every transcription
      step.  *breathing-only*: only find local breathing neighbors at every
      transcription step.""")
  parser.add_argument("--min-breathing", type=int, default=6, 
      metavar='<int>', 
      help="""Minimum number of freed bases during helix breathing.  Breathing
      helices can vary greatly in length, starting with at least two
      base-pairs. This parameter defines the minimum amount of bases freed by
      helix breathing. For example, 6 corresponds to a stack of two base-pairs
      and a loop region of 2 nucleotides. If less bases are freed and there
      exists a nested stacked helix, this helix is considered to breathe as
      well.""")

  # Advanced plotting parameters
  parser.add_argument("--t-lin", type=int, default=30, metavar='<int>',
      help="""Evenly space output *t-lin* times during transcription on a
      linear time scale.""")
  parser.add_argument("--t-log", type=int, default=300, metavar='<int>',
      help="""Evenly space output *t-log* times after transription on a
      logarithmic time scale.""")
  #NOTE: Needs to be adjusted after treekin dependency is gone...
  parser.add_argument("--t0", type=float, default=1e-4, metavar='<flt>',
      help=argparse.SUPPRESS)

  # More supported library parameters
  ril.argparse_add_arguments(parser, start=True, stop=True, 
      tmpdir=True, name=True, verbose=True)

  parser.add_argument("--k0", type=float, default=2e5, metavar='<flt>',
      help="""Arrhenius rate constant. Adjust the rate constant k0 of the the
      Arrhenius equation to match experimentally confirmed folding
      time-scales.""")

  # Plotting tools (DrForna, matplotlib, xmgrace)
  parser.add_argument("--drffile", action="store_true",
      help="Write DrForna output to a file: {--name}.drf") 
  parser.add_argument("--pyplot", action="store_true",
      help="""Plot the simulation using matplotlib. Interpret the legend using
      STDOUT or --logfile""")
  parser.add_argument("--xmgrace", action="store_true",
      help="""Plot the simulation for xmgrace visualization. Interpret the
      legend using STDOUT or --logfile""")

  # Logging and STDOUT 
  parser.add_argument("--logfile", action="store_true",
      help="Write verbose information to a file: {--name}.log") 
  parser.add_argument("--stdout", default='log', action = 'store',
      choices=('log', 'drf'),
      help="""Choose one of two STDOUT formats: *log*: prints {--verbose}
      logging information about the cotranscriptional folding progress, *drf*:
      prints the input format for DrForna to STDOUT.""")

  return

def main(args):
  """ DrTransformer - cotranscriptional folding """

  (name, fullseq) = ril.parse_vienna_stdin(sys.stdin)

  # Adjust arguments, prepare simulation
  if args.name == '' : 
    args.name = name
  else :
    name = args.name

  # Adjust filehandle-stuff
  _drffile = None
  if args.drffile :
    _drffile = name + '.drf'
    if args.stdout == 'drf' :
      args.stdout = ''
  elif args.stdout == 'drf' :
    _drffile = '-'

  _logfile = None
  if args.logfile:
    _logfile = name + '.log'
    if args.stdout == 'log' :
      args.stdout = ''
  elif args.stdout == 'log' :
    _logfile = '-'

  if args.tmpdir :
    _tmpdir = args.tmpdir 
    if not os.path.exists(_tmpdir):
      os.makedirs(_tmpdir)
  else :
    _tmpdir = tempfile.mkdtemp(prefix='DrTrafo')

  # Adjust simulation parameters
  args._RT=0.61632077549999997
  if args.temperature != 37.0 :
    kelvin = 273.15 + args.temperature
    args._RT = (args._RT/310.15)*kelvin

  if args.stop is None : 
    args.stop = len(fullseq)+1
  else :
    fullseq = fullseq[0:args.stop-1]

  if args.tX < args.t8 :
    raise ValueError('Simulation time after transcription --tX must be >= --t8')

  # Minrate specifies the lowest accepted rate for simulations (sec^-1)
  # it can be directly converted into a activation energy that results this rate
  args.maxdG = -args._RT * math.log(args.min_rate)
  #print args.min_rate, '=>', args.maxdG

  ############################
  # Start with DrTransformer #
  ############################

  # Set model details.
  vrna_md = RNA.md()
  vrna_md.noLP = 1
  vrna_md.temperature = args.temperature
  #vrna_md.circ = 2
  #vrna_md.dangles = 2

  #RNA.cvar.noLonelyPairs = 1
  #RNA.cvar.temperature = args.temperature

  vrna_fc = RNA.fold_compound(fullseq, vrna_md)
  #fc.contraints_add('.....(..)', RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP)

  check_version(args.treekin, ril._MIN_TREEKIN_VERSION)

  if _drffile :
    with smart_open(_drffile, 'w') as dfh :
      dfh.write("id time conc struct energy\n")

  if args.pyplot or args.xmgrace:
    all_courses = []

  if _logfile :
    with smart_open(_logfile, 'w') as lfh :
      lfh.write("# {} \n".format(fullseq))
      lfh.write("# ID, Structure, Energy, Occupancy\n")

  # initialize a directed conformation graph
  CG = nx.DiGraph(
      full_sequence = fullseq, 
      model_details = vrna_md,
      fold_compound = vrna_fc,
      transcript_length = None,
      total_time = 0,
      seqid = 0)
  saddles = dict()

  norm, plusI, expo = 0, 0, 0
  for tlen in range(args.start, args.stop) :
    CG.graph['transcript_length']=tlen
    seq = fullseq[0:tlen]

    # Expand Graph
    nn = expand_graph(CG, saddles, args, mode=args.structure_search_mode)
    #print """ {} new nodes """.format(nn), CG.graph['seqid'], "total nodes"

    if args.pyplot or args.xmgrace:
      ttt = CG.graph['total_time']
      if ttt == 0 :
        all_courses.extend([[] for i in range(nn)])
      else :
        all_courses.extend([[(ttt,0)] for i in range(nn)])
      # DO **NOT** DO IT THIS WAY: all_courses.extend([ [] ] * nn )

    # Simulate
    _fname = _tmpdir+'/'+name+'-'+str(tlen)
    _t0 = args.t0 if args.t0 > 0 else 1e-6
    _t8 = args.tX if tlen == args.stop-1 else args.t8
    (t_lin, t_log) = (None, args.t_log) if tlen == args.stop-1 else (args.t_lin, None)

    # produce input for treekin simulation
    with smart_open(_logfile, 'a') as lfh :
      [bfile, rfile, p0, nlist] = dump_conformation_graph(CG, seq, saddles, _fname, 
          logf=lfh, verb=(False or _logfile))

    dn,sr = 0,0
    if len(nlist) == 1 :
      # Fake Results for DrForna
      CG.graph['total_time'] += _t8
      if _drffile :
        with smart_open(_drffile, 'a') as dfh:
          ss = nlist[0][0]
          line = "{} {} {} {:s} {:6.2f}".format(CG.node[ss]['identity'], 
              CG.graph['total_time'], 1.0, ss[:len(seq)], CG.node[ss]['energy'])
          dfh.write(line + '\n')

      if args.pyplot or args.xmgrace:
        ss = nlist[0][0]
        ident = CG.node[ss]['identity']
        all_courses[ident].append((CG.graph['total_time'], 1.0))

    else :
      bfile = None # sometimes bfile causes a segfault, so let's leave it out.
      try: # - Simulate with treekin
        tfile, _ = ril.sys_treekin(_fname, seq, bfile, rfile, binrates=True,
            treekin=args.treekin, p0=p0, t0=_t0, ti=args.ti, t8=_t8, mpack=True,
            exponent=False, useplusI=False, force=True, verb=(args.verbose > 1))
        norm += 1
      except SubprocessError: 
        try : # - Simulate with treekin and --exponent
          tfile, _ = ril.sys_treekin(_fname, seq, bfile, rfile, binrates=True,
              treekin=args.treekin, p0=p0, t0=_t0, ti=args.ti, t8=_t8, 
              exponent=True, useplusI=False, force=True, verb=(args.verbose>0))
          expo += 1
        except SubprocessError: 
          try : # - Simulate with treekin and --useplusI
            tfile, _ = ril.sys_treekin(_fname, seq, bfile, rfile, binrates=True,
                treekin=args.treekin, p0=p0, t0=_t0, ti=args.ti, t8=_t8, 
                exponent=False, useplusI=True, force=True, verb=(args.verbose>0))
            plusI += 1
          except SubprocessError:
            if args.verbose > 1:
              print Warning("After {} nucleotides: treekin cannot find a solution!".format(tlen))
            # - Simulate with crnsimulator python package (slower)
            _odename = name+str(tlen)
            tfile = DiGraphSimulator(CG, _fname, nlist, p0, _t0, _t8, 
                t_lin = t_lin,
                t_log = t_log,
                jacobian=False, # faster!
                verb=(args.verbose>0))

      except ExecError, e:
        # NOTE: This is a hack to avoid treekin simulations in the first place
        _odename = name+str(tlen)
        tfile = DiGraphSimulator(CG, _fname, nlist, p0, _t0, _t8, 
            t_lin = t_lin,
            t_log = t_log,
            jacobian=False, # faster!
            verb=(args.verbose>1))

      # Get Results
      time_inc, iterations = get_stats_and_update_occupancy(
          CG, nlist, tfile)
      #print time_inc, iterations

      if args.pyplot or args.xmgrace or _drffile :
        for data in talk_generator(CG, nlist, tfile) :
          [id_, tt_, oc_, ss_, en_] = data
          if args.pyplot or args.xmgrace :
            all_courses[id_].append((tt_,oc_))
          if _drffile :
            with smart_open(_drffile, 'a') as dfh:
              dfh.write("{} {} {} {:s} {:6.2f}\n".format(*data))

      CG.graph['total_time'] += time_inc
      
      # Prune
      dn,sr = graph_pruning(CG, nlist, saddles, args)

    if args.verbose :
      print "# Transcripton length: {}. Initial graph size: {}. ".format(tlen, len(nlist)), 
      print "Deleted {} nodes, {} still reachable.".format(dn, sr)

  #if args.verbose >= 1:
  #  print "Treekin stats: {} default success, {} expo success, {} plusI success".format(
  #      norm, expo, plusI)

  if args.logfile or args.stdout == 'log':
    with smart_open(_logfile, 'a') as lfh :
      dump_conformation_graph(CG, CG.graph['full_sequence'], saddles, None,
          logf=lfh, verb=True)

  # CLEANUP the /tmp/directory
  if not args.tmpdir :
    shutil.rmtree(_tmpdir)

  # Plot results
  if args.pyplot:
    plot_simulation(all_courses, args)

  if args.xmgrace:
    plot_xmgrace(all_courses, args)

  return

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
      #formatter_class=argparse.RawTextHelpFormatter,
      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
      #formatter_class=argparse.MetavarTypeHelpFormatter,
      description='echo sequence | %(prog)s [options]')

  add_drtrafo_args(parser)

  args = parser.parse_args()

  main(args)

