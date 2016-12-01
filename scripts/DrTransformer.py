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
import collections as c
import tempfile
import shutil
import contextlib

import ribolands as ril
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
  """
  maxdG=args.maxdG
  k0 = args.k0
  fpath = args.findpath
  _RT = args._RT
  fullseq = CG.graph['full_sequence']


  saddleE = None
  if (s1, s2) in saddles :
    saddleE = saddles[(s1,s2)]
  else :
    saddleE = float(RNA.find_saddle(fullseq, s1, s2, fpath))/100

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
    k_12 = k0 * math.e**(-dG_1s/_RT) if dG_1s > 0 else k0
    k_21 = k0 * math.e**(-dG_2s/_RT) if dG_2s > 0 else k0

    CG.add_weighted_edges_from([(s1, s2, k_12)])
    CG.add_weighted_edges_from([(s2, s1, k_21)])

    #print "#Added Edge:", s1, s2, "({}, {:g}, {:g})".format(valid, k_12, k_21)

  return valid

def dump_conformation_graph(CG, seq, name, logf=sys.stdout,verb=False) :
  """ Make a barriers + rates output from the current conformation graph """

  sorted_nodes = sorted(CG.nodes(data=True), 
      key=lambda x: x[1]['energy'], reverse=False)

  def remove_inactive(snodes) :
    nnodes = []
    for (n,d) in snodes:
      if d['active'] == True :
        nnodes.append((n,d))
    return nnodes

  sorted_nodes = remove_inactive(sorted_nodes)

  if name :
    bfile = name+'.bar'
    rfile = name+'.rts'
    p0 = []
    with open(bfile, 'w') as bar, open(rfile, 'w') as rts :
      bar.write("     {}\n".format(seq))
      for e, (ni, data) in enumerate(sorted_nodes) :
        bar.write("{:4d} {} {:6.2f}\n".format(
          e+1, ni[:len(seq)], data['energy']))
        if verb :
          line = "{:4d} {:4d} {} {:6.2f} {:6.4f} (ID = {:d})\n".format(
              CG.graph['transcript_length'], e+1, ni[:len(seq)], 
              data['energy'], data['occupancy'], data['identity'])
          logf.write(line)

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
  else :
    line = "Distribution of structures at the end:\n"
    for e, (ni, data) in enumerate(sorted_nodes) :
      line += "LAST {:4d} {} {:6.2f} {:6.4f} (ID = {:d})\n".format(e+1, 
          ni[:len(seq)], data['energy'], data['occupancy'], data['identity'])
    logf.write(line)

    return 
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

def graph_pruning(CG, sorted_nodes, saddles, args) :
  """ Delete nodes or report them as still reachable. """
  cutoff = args.cutoff

  deleted_nodes = 0
  still_reachables = 0

  for ni, data in reversed(sorted_nodes) :
    #print ni, data
    if data['occupancy'] < cutoff :
      nbrs = sorted(CG.successors(ni), 
          key=lambda x: CG.node[x]['energy'], reverse=False)

      def remove_inactive(CG, snodes) :
        """ there must be a more beautiful way to do this ... """
        nnodes = []
        for n in snodes:
          if CG.node[n]['active'] == True :
            nnodes.append((n))
        return nnodes
      nbrs = remove_inactive(CG, nbrs)

      best, been = nbrs[0], CG.node[nbrs[0]]['energy']

      if been > data['energy'] :
        still_reachables += 1
        continue

      for e, nbr in enumerate(nbrs[1:]) :
        always_true = add_transition_edges(CG, saddles, args, nbr, best, ni)
        if always_true is False :
          sys.exit('over and out')

      if True: # hack here to keep all nodes
        CG.node[ni]['active']=False
        CG.node[ni]['occupancy']=0.0
        deleted_nodes += 1

  return deleted_nodes, still_reachables

def expand_graph(CG, saddles, args, mfe_only=False):
  """ Find new neighbors and add them to the Conformation Graph

  The function is devided into two parts. 1) The current mfe structure
  is connected to all present structures, 2) The conformation graph is
  expanded using helix-breathing.

  :param CG: Conformation Graph (NetworkX)
  :param saddles: dictionary of all previous findpath runs
  :param args: commandline arguments and other global variables
    (using: cutoff, verbose)
  :param mfe_only: only use current mfe as potential new neighbor

  :return: Number of new nodes

  """
  cutoff= args.cutoff
  verb  = args.verbose

  csid = CG.graph['seqid']
  fseq = CG.graph['full_sequence']
  tlen = CG.graph['transcript_length']
  seq = fseq[0:tlen]

  # Add MFE
  ss, mfe = RNA.fold(seq)
  future = '.' * (len(fseq)-tlen)
  ss = ss + future
  #print >> sys.stderr, "{}\n{} {:6.2f}".format(seq, ss, mfe)

  # If there is no node bec we are in the beginning, add the node,
  # otherwise, go through all nodes and try to add transition edges

  if nx.number_of_nodes(CG) == 0 :
    en = round(RNA.energy_of_structure(fseq, ss, 0), 2)
    CG.add_node(ss, energy=en, occupancy=1.0, 
        identity=CG.graph['seqid'], active=True)
    CG.graph['seqid'] += 1
  else :
    for ni in CG.nodes() :
      if CG.node[ni]['active'] == False : continue
      if ni == ss or CG.has_edge(ni,ss) : continue

      if CG.has_node(ss): # from a previous iteration
        if add_transition_edges(CG, saddles, args, ni, ss):
          CG.node[ss]['active'] = True # in case it was there but inactive
      elif add_transition_edges(CG, saddles, args, ni, ss) :
        en = round(RNA.energy_of_structure(fseq, ss, 0), 2)
        CG.node[ss]['active'] = True
        CG.node[ss]['energy'] = en
        CG.node[ss]['occupancy'] = 0.0
        CG.node[ss]['identity'] = CG.graph['seqid']
        CG.graph['seqid'] += 1

  if not CG.has_node(ss) or CG.node[ss]['active'] is False:
    print "# WARNING: ", ss, "[mfe secondary structure could not be connected]"

  if mfe_only is True :
    pass
  else :
    """ do the helix breathing graph expansion """
    for ni, data in CG.nodes_iter(data=True):
      if data['active'] == False : continue
      en  = data['energy']
      occ = data['occupancy']
      if occ < cutoff : continue

      ss = ni[0:len(seq)]

      opened = open_breathing_helices(seq, ss)
      for onbr in opened :
        nbr = fold_exterior_loop(seq, onbr)
        future = '.' * (len(ni) - len(nbr))
        nbr += future

        if ni == nbr or CG.has_edge(ni, nbr):
          continue

        if CG.has_node(nbr):
          if add_transition_edges(CG, saddles, args, ni, nbr):
            CG.node[nbr]['active'] = True # in case it was there but inactive
        elif add_transition_edges(CG, saddles, args, ni, nbr) :
          enbr = round(RNA.energy_of_structure(fseq, nbr, 0), 2)
          CG.node[nbr]['energy'] = enbr
          CG.node[nbr]['active'] = True
          CG.node[nbr]['occupancy'] = 0.0
          CG.node[nbr]['identity'] = CG.graph['seqid']
          CG.graph['seqid'] += 1
        else :
          """# WARNING: Could not add transition edge!"""
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
      else :
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

def fold_exterior_loop(seq, con, exterior_only=True):
  """ Constrained folding 
  
  The default behavior is "exterior_only", which replaces all constrained
  helices with short 'NNN' stretches at the sequence level. This reduces 
  the sequence length (n) and therefore the runtime O(n^3)
  
  :param seq: RNA sequence
  :param con: constraint
  :param exterior_only: only fold the extior loop region

  :return: secondary structure
  """

  if exterior_only :
    spacer = 'NNN'
    pt = ril.make_pair_table(con, base=0)
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

def talk_generator(CG, sorted_nodes, tfile, repl):
  """ Generator function that yields the time course information from the
  treekin output file.
  
  :param CG: Conformation Graph (networkx)
  :param sorted_nodes: a list of nodes sorted by their energy
  :param tfile: current treekin-output file
  :param repl: a factor to reduce the plot size, see --repl <float> option

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

          # can we actually skip the point?
          if repl and prevcourse and prevcourse[e+1] > 0 \
              and line is not tknlines[-2]:
            y1 = min(prevcourse[e+1], course[e+1])
            y2 = max(prevcourse[e+1], course[e+1])
            dy = math.log(y2) - math.log(y1)
            if dy < repl : continue

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
@line linewidth .1
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
  import matplotlib.pyplot as plt
  t8 = args.t8
  stop = args.stop
  cutoff = args.cutoff

  name  = args.name
  title = args.name

  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  #ax.set_xscale('log')
  for e, course in enumerate(all_in) :
    if course == [] : continue
    t, o = zip(*course)
    # determine which lines are part of the legend,
    # like this, it is only those that are populated 
    # at the end of transcription
    if t[-1] > t8*stop :
      #print e, t[-1], o[-1]
      p, = ax.plot(t, o, '--')
      p.set_label("ID {:d}".format(e))
    else :
      p, = ax.plot(t, o, '-')

  fig.set_size_inches(7,3)
  fig.text(0.5,0.95, title, ha='center', va='center')
  plt.xlabel('time [seconds]', fontsize=11)
  plt.ylabel('occupancy [mol/l]', fontsize=11)
  for tlen in range(args.start, args.stop) :
    plt.axvline(x=tlen*t8, linewidth=0.05, color='black', linestyle='-')

  # """ Add ticks for 1 minute, 1 hour, 1 day, 1 year """
  # plt.axvline(x=60, linewidth=1, color='black', linestyle='--')
  # plt.axvline(x=3600, linewidth=1, color='black', linestyle='--')
  # plt.axvline(x=86400, linewidth=1, color='black', linestyle='--')
  # plt.axvline(x=31536000, linewidth=1, color='black', linestyle='--')
  plt.legend()

  pfile = name+'.pdf'
  plt.savefig(pfile, bbox_inches='tight')

  return

def add_drtrafo_args(parser):
  """ A collection of arguments that are used by DrTransformer """
  parser.add_argument("--findpath", type = int, default = 10, metavar='<int>',
      help="Specify search width for *findpath* heuristic") 
  parser.add_argument("--minrate", type = float, default = 1e-10, metavar='<flt>',
      help="Specify minmum rate to accept a new neighbor")

  ril.argparse_add_arguments(parser, cutoff=True, start=True, stop=True, 
      tmpdir=True, name=True, verbose=True, k0=True, tX=True,
      temperature=True, treekin=True)

  parser.add_argument("--plot_cutoff", type=float, metavar='<flt>', default=0.02,
      help="Occupancy cutoff for final figure")
  #parser.add_argument("--plot_title", default='')
  #parser.add_argument("--plot_linlog", action="store_true",
  #    help="Divide x-axis into lin and log at transcription stop")

  
  parser.add_argument("--repl", type=float, default=None, metavar='<flt>',
      help="Logarithmically reduce output size")
  parser.add_argument("--logfile", action="store_true",
      help="Write verbose information into name.log") 
  parser.add_argument("--drffile", action="store_true",
      help="Write DrForna output into name.drf") 
  parser.add_argument("--pyplot", action="store_true",
      help="Plot the simulation using matplotlib. Interpret the legend \
          using the *log* output")
  parser.add_argument("--xmgrace", action="store_true",
      help="Print a plot for xmgrace. " + \
          "Interpret the legend using the *log* output")
  parser.add_argument("--stdout", default='log', action = 'store',
      metavar='<str>', help="Choose the stdout format: <log, drf>")

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

  # Minrate specifies the lowest accepted rate for simulations (sec^-1)
  # it can be directly converted into a activation energy that results this rate
  args.maxdG = -args._RT * math.log(args.minrate)
  #print args.minrate, '=>', maxdG

  ############################
  # Start with DrTransformer #
  ############################
  RNA.cvar.noLonelyPairs = 1
  RNA.cvar.temperature = args.temperature

  if _drffile :
    with smart_open(_drffile, 'w') as dfh :
      dfh.write("id time conc struct energy\n")

  if args.pyplot or args.xmgrace:
    all_courses = []

  if _logfile :
    with smart_open(_logfile, 'w') as lfh :
      lfh.write(fullseq + "\n")
      lfh.write("# ID, Structure, Energy, Occupancy\n")

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
    nn = expand_graph(CG, saddles, args, mfe_only=False)
    #print """ {} new nodes """.format(nn), CG.graph['seqid'], "total nodes"

    if args.pyplot or args.xmgrace:
      all_courses.extend([[] for i in range(nn)])
      # Just so that I will remember... 
      # DO **NOT** DO IT THIS WAY: all_courses.extend([ [] ] * nn )
      # and if you are a python geek reading this, please tell me when or why
      # you would ever allocate space like that!

    # Simulate
    _fname = _tmpdir+'/'+name+'-'+str(tlen)
    _t8 = args.tX if tlen == args.stop-1 else args.t8

    # produce input for treekin simulation
    with smart_open(_logfile, 'a') as lfh :
      [bfile, rfile, p0, nlist] = dump_conformation_graph(CG, seq, _fname, 
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
      # - Simulate with treekin
      try:
        tfile, _ = ril.sys_treekin(_fname, seq, bfile, rfile, 
            treekin=args.treekin, p0=p0, t0=args.t0, ti=args.ti, t8=_t8, 
            useplusI=True, force=True, verb=(args.verbose > 1))
      except RuntimeError:
        #try :
        #  tfile = ril.sys_treekin(_fname, seq, bfile, rfile, 
        #      treekin=args.treekin, p0=p0, t0=args.t0, ti=args.ti, t8=_t8, 
        #      useplusI=False, force=True, verb=True)
        #except RuntimeError:
          print "Abort after", tlen, "nucleotides:", \
              "treekin cannot find a solution, sorry"
          raise SystemExit

      # Get Results
      time_inc, iterations = get_stats_and_update_occupancy(CG, nlist, tfile)
      #print time_inc, iterations

      if args.pyplot or args.xmgrace or _drffile :
        for data in talk_generator(CG, nlist, tfile, args.repl) :
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
      print "# Deleted {} nodes, {} still reachable.".format(dn, sr)

  if args.logfile or args.stdout == 'log':
    with smart_open(_logfile, 'a') as lfh :
      dump_conformation_graph(CG, CG.graph['full_sequence'], None,
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

