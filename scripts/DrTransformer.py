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

import RNA
import ribolands as ril
from ribolands.syswraps import SubprocessError, ExecError, check_version, VersionError
from ribolands.crnwrapper import DiGraphSimulator
from ribolands.trafo import ConformationGraph

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

  parser.add_argument("--mpack-method", metavar='<str>', action='store',
      choices=(['FLOAT128', 'LD', 'QD', 'DD']),
      help="""Increase the precision of treekin simulations. Requires a development
      version of treekin with mpack support.""")

  # Common parameters
  parser.add_argument("--ti", type=float, default=1.2, metavar='<flt>',
      help="""Output-time increment of treekin solver (t1 * ti = t2).""")
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

  parser.add_argument("--minh", type=float, default=0.01, metavar='<flt>',
      help="Merge structures separated by a barrier smaller than minh.")

  parser.add_argument("--soft-minh", type=float, default=0, metavar='<flt>',
      help="""Merge structures separated by a barrier smaller than minh *for
      visualzation only*. The dynamics will still be caculated based on the
      more detailed network. This parameter will only have an effect if it is
      higher than --minh.""")

  parser.add_argument("--min-rate", type = float, default = None, metavar='<flt>',
      help="""Minmum rate to accept a new structure as neighboring
      conformation.""")

  parser.add_argument("--maxh", type=float, default=None, metavar='<flt>',
      help=argparse.SUPPRESS)
      #help="""Do *not* merge structures separated by a barrier higher than maxh
      #during graph pruning!""")

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
  _RT=0.61632077549999997
  if args.temperature != 37.0 :
    kelvin = 273.15 + args.temperature
    _RT = (_RT/310.15)*kelvin

  if args.stop is None : 
    args.stop = len(fullseq)+1
  else :
    fullseq = fullseq[0:args.stop-1]

  if args.tX < args.t8 :
    raise ValueError('Simulation time after transcription --tX must be >= --t8')

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # Adjust the simulation window for treekin:
  #
  #   k0           
  #   2e5  /s       50 nuc/s               1e-inf /s
  #   4000 /t8      1  nuc/t8              1e-inf /t8
  #   |------$------|-------------------$----------->
  #          k-fast                     k-slow
  #          --minh                     --maxh
  #          |------|---------------|--->
  #          t0     t8 0.02 s       tX = 86400 s
  #   <----->    simulation             <---------->
  #   instant                             rate = 0
  #
  # (1) The rate of a spontaneous folding event 
  # (>= k-fast) has to be much faster than the 
  # rate of chain elongation (kx).
  #
  # (2) The rate of an effectively 0 folding 
  # event (< k-slow) has to be much slower than 
  # the rate of chain elongation (kx), it 
  # should also be much slower than the 
  # post-transcriptional simulation time --tX
  #
  # Parameters:
  # k0 = maximum folding rate /s
  # t8 = time for chain elongation 
  # tX = post-transcriptional simulation time
  # t0 = first simulation output time (<< t8)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  if args.min_rate is not None:
    args.maxdG = -_RT * math.log(args.min_rate/args.k0)
  else :
    args.maxdG = 0 # off

  mhrate = args.k0 * math.exp(-args.minh/_RT)

  # Compare --minh, --min-rate, --t8 and --tX :
  if args.verbose:
    print '# {} dG minh => {} /s rate and {} s time at k0 = {}'.format(
        args.minh, mhrate, 1/mhrate, args.k0)

    if args.min_rate is not None:
      print '# {} /s min-rate => {} dG barrier and {:g} s time at k0 = {}'.format(
          args.min_rate, args.maxdG, 1/args.min_rate, args.k0)

    # Minrate specifies the lowest accepted rate for simulations (sec^-1)
    # it can be directly converted into a activation energy that results this rate
    print '# Settings correspond to potential new paramters:'
    print '# --t-fast: {} s, --t-ext {} s, --t-end {} s, --t-slow {} s'.format(
        1/mhrate, args.t8, args.tX, 1/args.min_rate if args.min_rate is not None else 'inf')

  if 1/mhrate * 10 > args.t8 :
    raise Exception('Conflicting Settings: rate for equilibration must be much faster than for nucleotide extension.')
    
  if args.min_rate is not None and args.tX * 100 > 1/args.min_rate :
    raise Exception('Conflicting Settings: 1/--min-rate should be much longer than the simulation time --tX.')


  check_version(args.treekin, ril._MIN_TREEKIN_VERSION)

  def versiontuple(rv):
    return tuple(map(int, (rv.split("."))))
  if versiontuple(RNA.__version__) < versiontuple(ril._MIN_VIENNARNA_VERSION):
    raise VersionError('ViennaRNA', RNA.__version__, ril._MIN_VIENNARNA_VERSION)

  ############################
  # Start with DrTransformer #
  ############################

  # Set model details.
  vrna_md = RNA.md()
  vrna_md.noLP = 1
  vrna_md.temperature = args.temperature

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
  CG = ConformationGraph(fullseq, vrna_md)
  CG.occupancy_cutoff = args.occupancy_cutoff
  CG.findpath_width = args.findpath_width
  CG._maxdG = args.maxdG
  CG._k0 = args.k0

  CG._transcript_length = args.start

  from datetime import datetime
  mytime = datetime.now()

  with smart_open(_logfile, 'w') as lfh :
    CG.logfile = lfh
    norm, plusI, expo, fail = 0, 0, 0, 0
    for tlen in range(args.start, args.stop) :
      nn = CG.expand(exp_mode=args.structure_search_mode)
      #print """ {} new nodes """.format(nn), CG._nodeid, "total nodes"

      mn = CG.coarse_grain(minh = args.minh)
      if args.verbose :
        print "# Merged {} nodes after expanding {} new nodes.".format(len(mn), nn)

      if args.pyplot or args.xmgrace:
        ttt = CG._total_time
        if ttt == 0 :
          all_courses.extend([[] for i in range(nn)])
        else :
          all_courses.extend([[(ttt,0)] for i in range(nn)])
        # DO **NOT** DO IT THIS WAY: all_courses.extend([ [] ] * nn )

      ############
      # Simulate #
      ############
      _fname = _tmpdir+'/'+name+'-'+str(tlen)
      _t0 = args.t0 if args.t0 > 0 else 1e-6
      _t8 = args.tX if tlen == args.stop-1 else args.t8
      (t_lin, t_log) = (None, args.t_log) if tlen == args.stop-1 else (args.t_lin, None)

      # produce input for treekin simulation
      [bfile, rfile, p0, nlist] = CG.get_simulation_files_tkn(_fname)

      for e, (ni, data) in enumerate(nlist, 1) :
          lfh.write("{:4d} {:4d} {} {:6.2f} {:6.4f} (ID = {:d})\n".format(
              tlen, e, ni[:tlen], data['energy'], data['occupancy'], data['identity']))

      dn, sr, rj = 0, 0, 0
      if len(nlist) == 1 :
        # Fake simulation results for DrForna
        CG._total_time += _t8
        if _drffile :
          with smart_open(_drffile, 'a') as dfh:
            ss = nlist[0][0]
            line = "{} {} {} {:s} {:6.2f}".format(CG.node[ss]['identity'], 
                CG._total_time, 1.0, ss[:len(seq)], CG.node[ss]['energy'])
            dfh.write(line + '\n')

        if args.pyplot or args.xmgrace:
          ss = nlist[0][0]
          ident = CG.node[ss]['identity']
          all_courses[ident].append((CG._total_time, 1.0))

      else :
        seq = ''
        bfile = None # sometimes bfile causes a segfault, so let's leave it out.
        try: # - Simulate with treekin
          tfile, _ = ril.sys_treekin(_fname, seq, bfile, rfile, binrates=True,
              treekin=args.treekin, p0=p0, t0=_t0, ti=args.ti, t8=_t8, 
              mpack_method=args.mpack_method,
              exponent=False, useplusI=False, force=True, verb=(args.verbose > 1))
          norm += 1
        except SubprocessError: 
          try : # - Simulate with treekin and --exponent
            tfile, _ = ril.sys_treekin(_fname, seq, bfile, rfile, binrates=True,
                mpack_method=args.mpack_method,
                treekin=args.treekin, p0=p0, t0=_t0, ti=args.ti, t8=_t8, 
                exponent=True, useplusI=False, force=True, verb=(args.verbose>0))
            expo += 1
          except SubprocessError: 
            try : # - Simulate with treekin and --useplusI
              tfile, _ = ril.sys_treekin(_fname, seq, bfile, rfile, binrates=True,
                  mpack_method=args.mpack_method,
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
              fail += 1

        except ExecError, e:
          # NOTE: This is a hack to avoid treekin simulations in the first place
          _odename = name+str(tlen)
          tfile = DiGraphSimulator(CG, _fname, nlist, p0, _t0, _t8, 
              t_lin = t_lin,
              t_log = t_log,
              jacobian=False, # faster!
              verb=(args.verbose>1))

        # Get Results
        time_inc, iterations = CG.update_occupancies_tkn(tfile, nlist)
        #print time_inc, iterations

        softmap = dict()
        if args.soft_minh and args.soft_minh > args.minh :
          copyCG = CG.graph_copy()
          softmap = copyCG.coarse_grain(minh=args.soft_minh)
          del copyCG

        if args.pyplot or args.xmgrace or _drffile :
          for data in CG.sorted_trajectories_iter(nlist, tfile, softmap) :
            [id_, tt_, oc_, ss_, en_] = data
            if args.pyplot or args.xmgrace :
              all_courses[id_].append((tt_,oc_))
            if _drffile :
              with smart_open(_drffile, 'a') as dfh:
                dfh.write("{} {} {} {:s} {:6.2f}\n".format(*data))

        CG._total_time += time_inc
        
        # Prune
        dn, sr, rj = CG.prune(maxh=args.maxh, mocca=2)

      if args.verbose :
        print "# Transcripton length: {}. Initial graph size: {}. ".format(tlen, len(nlist)) 
        print "#  Deleted {} nodes, {} still reachable, {} rejected deletions.".format(dn, sr, rj)
        print "#  Computation time: {}".format(datetime.now()-mytime)
        print tlen, (datetime.now()-mytime).total_seconds(), len(nlist), dn, sr, rj, ril.trafo.PROFILE['findpath-calls'],
        print ril.trafo.PROFILE['mfe'], ril.trafo.PROFILE['hb'], ril.trafo.PROFILE['feature'], ril.trafo.PROFILE['cogr'], ril.trafo.PROFILE['prune']

        mytime = datetime.now()
        ril.trafo.PROFILE['findpath-calls'] = 0
        ril.trafo.PROFILE['mfe'] = 0
        ril.trafo.PROFILE['hb'] = 0
        ril.trafo.PROFILE['feature'] = 0
        ril.trafo.PROFILE['cogr'] = 0
        ril.trafo.PROFILE['prune'] = 0

    if args.verbose >= 1:
      print "# Treekin stats: {} default success, {} expo success, {} plusI success, {} fail".format(
          norm, expo, plusI, fail)

    lfh.write("Distribution of structures at the end:\n")
    lfh.write("          {}\n".format(CG.transcript))
    for e, (ni, data) in enumerate(CG.sorted_nodes(), 1) :
        lfh.write("{:4d} {:4d} {} {:6.2f} {:6.4f} (ID = {:d})\n".format(
            tlen, e, ni[:tlen], data['energy'], data['occupancy'], data['identity']))

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

