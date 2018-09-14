#
# DrTransformer.py -- co-transcriptional folding.
#
# written by Stefan Badelt (stef@tbi.univie.ac.at)
#

from __future__ import absolute_import, division, print_function #, unicode_literals

from builtins import map
from builtins import zip
from builtins import str
from builtins import range

import os
import sys
import math
import argparse
import subprocess as s
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import shutil
import tempfile
from datetime import datetime

import RNA
import ribolands as ril
from ribolands.syswraps import SubprocessError, ExecError, check_version, VersionError
from ribolands.crnwrapper import DiGraphSimulator
from ribolands.trafo import TrafoLandscape

def plot_xmgrace(trajectories, filename):
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
    with open(filename, 'w') as gfh:
        gfh.write(head)
        for e, course in enumerate(trajectories):
            t, o = list(zip(*course))
            for i in range(len(t)):
                gfh.write("{:f} {:f}\n".format(t[i], o[i]))
            gfh.write("&\n")
    return


def plot_simulation(trajectories, 
        steps, 
        t_ext, 
        t_end,
        fpath,
        formats = ['pdf'],
        title = ''):
    """
    """

    # # Get the relevant arguments from args
    lin_time = steps * float(t_ext)
    log_time = lin_time + t_end if (lin_time + t_end) >= lin_time * 10 else lin_time * 10

    # Do the plotting
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylim([0, 1.01])
    ax.set_xscale('linear')

    # Make the second part of the plot logarithmic
    ax.set_xlim((0, lin_time))
    divider = make_axes_locatable(ax)
    axLog = divider.append_axes("right", size=2.5, pad=0, sharey=ax)
    axLog.set_xscale('log')
    axLog.set_xlim((lin_time + 0.00001, log_time))
    axLog.set_ylim([0, 1.01])
    axLog.yaxis.set_visible(False)

    for e, course in enumerate(trajectories):
        if course == []:
            continue
        t, o = list(zip(*course))
        # Determine which lines are part of the legend:
        # like this, it is only those that are populated
        # at the end of transcription and if they reach
        # an occupancy of 10% or higher
        if t[-1] > lin_time:
            p, = ax.plot(t, o, '-', lw=1.5)
            L, = axLog.plot(t, o, '-', lw=1.5)
            if max(o) >= 0.1:
                L.set_label("ID {:d}".format(e))
        else:
            p, = ax.plot(t, o, '-', lw=0.5)
            L, = axLog.plot(t, o, '-', lw=0.5)

    fig.set_size_inches(7, 3)
    fig.text(0.5, 0.95, title, ha='center', va='center')

    # for tlen in range(args.stop-args.start) :
    #  ax.axvline(x=tlen*t8, linewidth=0.01, color='black', linestyle='--')

    # Add ticks for 1 minute, 1 hour, 1 day, 1 year
    axLog.axvline(x=60, linewidth=1, color='black', linestyle='--')
    axLog.axvline(x=3600, linewidth=1, color='black', linestyle='--')
    axLog.axvline(x=86400, linewidth=1, color='black', linestyle='--')
    axLog.axvline(x=31536000, linewidth=1, color='black', linestyle='--')
    plt.legend()

    ax.set_ylabel('occupancy [mol/l]', fontsize=11)
    ax.set_xlabel('time [s]', ha='center', va='center', fontsize=11)
    ax.xaxis.set_label_coords(.9, -0.15)

    for ending in formats:
        pfile = fpath + '.' + ending
        plt.savefig(pfile, bbox_inches='tight')

    return


def add_drtrafo_args(parser):
    """ A collection of arguments that are used by DrTransformer """

    environ = parser.add_argument_group('DrTransformer dependencies')
    output = parser.add_argument_group('DrTransformer output')
    trans = parser.add_argument_group('Transcription parameters')
    algo  = parser.add_argument_group('DrTransformer algorithm')

    ###################
    # Default options #
    ###################
    parser.add_argument('--version', action='version', version='%(prog)s ' + ril.__version__)
    parser.add_argument("-v", "--verbose", action='count', default=0,
            help="""Track process by writing verbose output to STDOUT during
            calculations. Use --logfile if you want to see *just* verbose
            information via STDOUT.""")

    ##############################
    # DrTransformer dependencies #
    ##############################
    environ.add_argument("--treekin", default='treekin', action='store', metavar='<str>', 
            help="Path to the *treekin* executable.")

    environ.add_argument("--mpack-method", action='store',
            choices=(['FLOAT128', 'LD', 'QD', 'DD']),
            help="""Increase the precision of treekin simulations. Requires a development
                    version of treekin with mpack support.""")

    ########################
    # DrTransformer output #
    ########################
    output.add_argument("--stdout", default=None, action='store',
            choices=('log', 'drf', 'OFF'),
            help="""Choose STDOUT formats to follow the cotranscriptional
            folding progress in real time: *log*: a human readable output
            format.  *drf*: DrForna visalization input format. *OFF*: actively
            suppress output. The default (None) switches between *OFF*, if --logfile is
            specified, or *log* otherwise.""")

    output.add_argument("--tmpdir", default='', action='store', metavar='<str>',
            help="""Specify path for storing temporary output files. These
            files will not be removed when the program terminates. Defaults to
            the default directory for temporary files on your System.""")

    output.add_argument("--outdir", default='', action='store', metavar='<str>',
            help="""Place regular output files, into this directory. Creates
            the directory if it does not exist. """)

    output.add_argument("--name", default='', metavar='<str>',
            help="""Name your output files, name the header of your plots, etc.
            this option overwrites the fasta-header.""")

    output.add_argument("--logfile", action="store_true",
            help="""Write verbose information to a file:
            {--outdir}/{--name}.log""")

    output.add_argument("--visualize", nargs='+', default='', action='store',
            choices=('pdf', 'svg', 'png', 'gr', 'drf'),
            help="""Plot the simulation using matplotlib (pdf, svg, png) and/or
            write an input file for xmgrace (gr) and/or write an input file for
            DrForna (drf). Interpret the legend using STDOUT or --logfile.
            Files: {--outdir}/{--name}.{--visualize}""")

    output.add_argument("--tinc", type=float, default=None, metavar='<flt>',
            #NOTE: Deprecated: use --t-inc!
            help=argparse.SUPPRESS)

    output.add_argument("--t-inc", type=float, default=1.2, metavar='<flt>',
            help="""Adjust the plotting time resolution via the time-increment of
            the solver (t1 * t-inc = t2).""")

    output.add_argument("--soft-minh", type=float, default=0, metavar='<flt>',
            #NOTE: Deprecated: use --v-minh!
            help=argparse.SUPPRESS)

    output.add_argument("--plot-minh", type=float, default=0, metavar='<flt>',
            help="""Reduce the resolution of visualized structures. This
            parameter merges structures which are separated by a barrier
            smaller than --plot-minh *for visualzation only*. The dynamics will
            still be caculated based on the more detailed network. This
            parameter can only be effectve if it leads to a stronger
            coarse-graining than the simulation paramter --t-fast.  
            
            This parameter corresponds to folding time as: 
            t=1/(k0*exp(-dG/RT)), where dG is the parameter --plot-minh.
            E.g.  --plot-minh: 5.1 kcal/mol ~= 0.02 s [default t-ext]. """)

    output.add_argument("--draw-graphs", action="store_true",
            #help="""Export every landscape as json file. Uses --tempdir. """)
            help=argparse.SUPPRESS)

    output.add_argument("--t-lin", type=int, default=30, metavar='<int>',
            #help="""Evenly space output *t-lin* times during transcription on a linear time scale.""")
            help=argparse.SUPPRESS)
    output.add_argument("--t-log", type=int, default=300, metavar='<int>',
            #help="""Evenly space output *t-log* times after transcription on a logarithmic time scale.""")
            help=argparse.SUPPRESS)

    output.add_argument("--t0", type=float, default=0, metavar='<flt>',
            help=argparse.SUPPRESS)

    ############################
    # Transcription parameters #
    ############################
    trans.add_argument("--t-ext", type=float, default=0.02, metavar='<flt>',
            help="""Transcription speed, i.e. time per nucleotide extension
            [seconds per nucleotide].""")
    trans.add_argument("--t-end", type=float, default=60, metavar='<flt>',
            help="Post-transcriptional simulation time [seconds].")
    trans.add_argument("-T", "--temperature", type=float, default=37.0, metavar='<flt>', 
            help="The temperature for ViennaRNA computations.")
    trans.add_argument("--start", type=int, default=1, metavar='<int>',
            help="Start transcription at this nucleotide.")
    trans.add_argument("--stop", type=int, default=None, metavar='<int>',
            help="Stop transcription before this nucleotide")

    ###########################
    # DrTransformer algorithm #
    ###########################
    algo.add_argument("--min-occupancy", type=float, default=0.01, metavar='<flt>',
            help="""Occupancy threshold to determine which structures are
            considered (relevant) when transcribing a new nucleotide.
            Probability threshold for secondary structure graph expansion.""")

    algo.add_argument("--t-fast", type=float, default=5e-6, metavar='<flt>',
            help="""Folding times faster than --t-fast are considered
            instantaneous.  Structural transitions that are faster than
            --t-fast are considerd part of the same macrostate. Directly
            translates to an energy barrier separating conforations using:
            dG=-RT*ln(1/t_fast/k0).""")

    algo.add_argument("--t-slow", type=float, default=None, metavar='<flt>',
            help="""Only accept new structures as neighboring conformations if
            the transition is faster than --t-slow. This parameter may be
            useful to prevent numeric instabilities, otherwise, better avoid
            it.""")

    algo.add_argument("--findpath-search-width", type=int, default=20, metavar='<int>',
            help="""Search width for the *findpath* heuristic. Higher values
            increase the chances to find energetically better transition state
            energies.""")

    algo.add_argument('--structure-search-mode', default='default',
            choices=('default', 'mfe-only', 'breathing-only'),
            help="""Specify one of three modes: *default*: find new secondary
            structures using both the current MFE structure and breathing
            neighbors.  *mfe-only*: only find the current MFE structure at
            every transcription step.  *breathing-only*: only find local
            breathing neighbors at every transcription step.""")

    algo.add_argument("--min-breathing", type=int, default=6, metavar='<int>',
            help="""Minimum number of freed bases during helix breathing.
            Breathing helices can vary greatly in length, starting with at
            least two base-pairs. This parameter defines the minimum amount of
            bases freed by helix breathing. For example, 6 corresponds to a
            stack of two base-pairs and a loop region of 2 nucleotides. If less
            bases are freed and there exists a nested stacked helix, this helix
            is considered to breathe as well.""")

    algo.add_argument("--detailed-pruning", action="store_true",
            help="""Calculate new direct path barriers during the pruning step.""")

    algo.add_argument("--k0", type=float, default=2e5, metavar='<flt>',
            help="""Arrhenius rate constant. Adjust the rate constant k0 of the
            Arrhenius equation to relate free energy changes to experimentally
            confirmed folding time in seconds.""")

    algo.add_argument("--mocca", type=int, default=None, metavar='<int>',
            # I actually forgot what this was about, but it's an over-the-thumb solution,
            # which was worth experimenting with during profiling.
            help=argparse.SUPPRESS)

    return

def write_output(data, stdout = False, fh = None):
    # Helper function to print data to filehandle, STDOUT, or both.
    if stdout:
        sys.stdout.write(data)
    if fh:
        fh.write(data)


def main(args):
    """ DrTransformer - cotranscriptional folding """
    (name, fullseq) = ril.parse_vienna_stdin(sys.stdin)

    # Argument deprecation Warnings
    if args.tinc:
        print('DEPRECATION WARNING: Argument --tinc is deprecated, please use --t-inc.')
        args.t_inc = args.tinc

    if args.soft_minh:
        print('DEPRECATION WARNING: Argument --soft-minh is deprecated, please use --v-minh.')
        args.plot_minh = args.soft_minh


    # Adjust arguments, prepare simulation
    if args.name == '':
        args.name = name
    else:
        name = args.name

    if os.path.split(args.name)[0]:
        raise NotImplementedError('Name must not contain filepath.')

    if args.outdir:
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        filepath = args.outdir + '/' + args.name
    else :
        filepath = args.name

    dfh = open(filepath + '.drf', 'w') if 'drf' in args.visualize else None
    lfh = open(filepath + '.log', 'w') if args.logfile else None

    args.pyplot = filter(lambda x: x != 'drf', args.visualize)

    # Adjust filehandle-stuff
    if args.stdout is None and lfh is None:
        args.stdout = 'log'

    if args.tmpdir:
        _tmpdir = args.tmpdir
        if not os.path.exists(_tmpdir):
            os.makedirs(_tmpdir)
    else:
        _tmpdir = tempfile.mkdtemp(prefix='DrTrafo')

    # Adjust simulation parameters
    _RT = 0.61632077549999997
    if args.temperature != 37.0:
        kelvin = 273.15 + args.temperature
        _RT = (_RT / 310.15) * kelvin

    if args.stop is None:
        args.stop = len(fullseq) + 1

    if args.t_end < args.t_ext:
        raise ValueError(
            'Simulation time after transcription --t-end must be >= --t-ext')

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

    # Compare --minh, --min-rate, --t8 and --tX :
    if args.verbose:
        dG_min = -_RT * math.log(1 / args.t_fast / args.k0)
        print('# --t-fast: {} s => {} kcal/mol barrier height and {} /s rate at k0 = {}'.format(
            args.t_fast, dG_min, 1/args.t_fast, args.k0))
        if args.t_slow is not None:
            dG_max = -_RT * math.log(1 / args.t_slow / args.k0)
            print('# --t-slow {} s => {} kcal/mol barrier height and {} /s rate at k0 = {}'.format(
                args.t_slow, dG_max, 1/args.t_slow, args.k0))

    if args.t_fast * 10 > args.t_ext:
        raise Exception("""Conflicting Settings: rate for equilibration must be
                much faster than for nucleotide extension.""")

    if args.t_slow is not None and args.t_end * 100 > args.t_slow:
        raise Exception("""Conflicting Settings: 1/--min-rate should be much
                longer than the simulation time --t_end.""")

    try:
        check_version(args.treekin, ril._MIN_TREEKIN_VERSION)
    except ExecError:
        pass

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

    if args.pyplot:
        all_courses = []

    # Write logging output
    if args.stdout == 'log' or lfh:
        fdata  = "# File generated using DrTransformer v{}\n".format(ril.__version__)
        fdata += "#\n"
        fdata += "# >{}\n# {} \n".format(name, fullseq)
        fdata += "#\n"
        fdata += "# Co-transcriptional folding paramters:\n"
        fdata += "# --t-ext: {} sec\n".format(args.t_ext)
        fdata += "# --t-end: {} sec\n".format(args.t_end)
        fdata += "# --temperature: {} C\n".format(args.temperature)
        fdata += "# --start: {}\n".format(args.start)
        fdata += "# --stop: {}\n".format(args.stop)
        fdata += "#\n"
        fdata += "# Algorithm paramters:\n"
        fdata += "# --min-occupancy: {}\n".format(args.min_occupancy)
        fdata += "# --t-fast: {} sec\n".format(args.t_fast)
        fdata += "# --t-slow: {} sec\n".format(args.t_slow)
        fdata += "# --findpath-search-width: {}\n".format(args.findpath_search_width)
        fdata += "# --structure-search-mode: {}\n".format(args.structure_search_mode)
        fdata += "# --min-breathing: {} nuc\n".format(args.min_breathing)
        fdata += "# --detailed-pruning: {}\n".format(args.detailed_pruning)
        fdata += "# --k0: {}\n".format(args.k0)
        fdata += "#\n"
        fdata += "#\n"
        fdata += "# Results:\n"
        fdata += "# Tanscription Step | Energy-sorted structure index | Structure | Energy | Occupancy | Structure ID (-> Plotting ID)\n"
        write_output(fdata, stdout=(args.stdout == 'log'), fh = lfh)

    # Write DrForna output
    if args.stdout == 'drf' or dfh:
        fdata = "id time conc struct energy\n"
        write_output(fdata, stdout=(args.stdout == 'drf'), fh = dfh)

    # initialize a directed conformation graph
    CG = TrafoLandscape(fullseq, vrna_md)
    CG.p_min = args.min_occupancy
    CG._k0 = args.k0
    CG.t_fast = args.t_fast
    CG.t_slow = args.t_slow

    CG.findpath_search_width = args.findpath_search_width

    CG._transcript_length = args.start-1

    if args.verbose:
        atime = datetime.now()

    # now lets start...
    norm, plusI, expo, fail = 0, 0, 0, 0
    for tlen in range(args.start, args.stop):
        nn = CG.expand(exp_mode=args.structure_search_mode)
        # print """ {} new nodes """.format(nn), CG._nodeid, "total nodes"

        mn = CG.coarse_grain()
        if args.verbose:
            print("# Merged {} nodes after expanding {} new nodes.".format(len(mn), nn))

        if args.pyplot:
            ttt = CG._total_time
            if ttt == 0:
                all_courses.extend([[] for i in range(nn)])
            else:
                all_courses.extend([[(ttt, 0)] for i in range(nn)])

        ############
        # Simulate #
        ############
        _fname = _tmpdir + '/' + name + '-' + str(tlen)
        _t0 = args.t0  # if args.t0 > 0 else 1e-6
        _t8 = args.t_end if tlen == args.stop - 1 else args.t_ext
        (t_lin, t_log) = (None, args.t_log) if tlen == args.stop - \
            1 else (args.t_lin, None)

        # produce input for treekin simulation
        [bfile, rfile, p0, nlist] = CG.get_simulation_files_tkn(_fname)

        softmap = dict()
        if args.plot_minh and args.plot_minh > CG._dG_min:
            copyCG = CG.graph_copy()
            softmap = copyCG.coarse_grain(dG_min=args.plot_minh)
            del copyCG

        if args.stdout == 'log' or lfh:
            for e, (ni, data) in enumerate(nlist, 1):
                sm = '-> {}'.format(CG.node[softmap[ni]]['identity']) if ni in softmap else ''
                fdata = "{:4d} {:4d} {} {:6.2f} {:6.4f} (ID = {:d} {:s})\n".format(
                    tlen, e, ni[:tlen], data['energy'], data['occupancy'], data['identity'], sm)
                write_output(fdata, stdout=(args.stdout == 'log'), fh = lfh)

        dn, sr, rj = 0, 0, 0
        if len(nlist) == 1:
            # Fake simulation results for DrForna
            CG._total_time += _t8
            if args.stdout == 'drf' or dfh:
                ss = nlist[0][0]
                fdata = "{} {} {} {:s} {:6.2f}\n".format(CG.node[ss]['identity'],
                        CG._total_time, 1.0, CG.transcript, CG.node[ss]['energy'])
                write_output(fdata, stdout=(args.stdout == 'drf'), fh = dfh)

            if args.pyplot:
                ss = nlist[0][0]
                ident = CG.node[ss]['identity']
                all_courses[ident].append((CG._total_time, 1.0))

        else:
            seq = ''
            # sometimes bfile causes a segfault, so let's leave it out.
            bfile = None
            try:  # - Simulate with treekin
                tfile, _ = ril.sys_treekin(_fname, seq, bfile, rfile, binrates=True,
                        treekin=args.treekin, p0=p0, t0=_t0, ti=args.t_inc, t8=_t8,
                        mpack_method=args.mpack_method,
                        exponent=False, useplusI=False, force=True, verb=(args.verbose > 1))
                norm += 1
            except SubprocessError:
                try:  # - Simulate with treekin and --exponent
                    tfile, _ = ril.sys_treekin(_fname, seq, bfile, rfile, binrates=True,
                            mpack_method=args.mpack_method,
                            treekin=args.treekin, p0=p0, t0=_t0, ti=args.t_inc, t8=_t8,
                            exponent=True, useplusI=False, force=True, verb=(args.verbose > 0))
                    expo += 1
                except SubprocessError:
                    try:  # - Simulate with treekin and --useplusI
                        tfile, _ = ril.sys_treekin(_fname, seq, bfile, rfile, binrates=True,
                                mpack_method=args.mpack_method,
                                treekin=args.treekin, p0=p0, t0=_t0,
                                ti=args.t_inc, t8=_t8, exponent=False,
                                useplusI=True, force=True,
                                verb=(args.verbose > 0))
                        plusI += 1
                    except SubprocessError:
                        if args.verbose > 1:
                            print("After {} nucleotides: treekin cannot find a solution!".format(tlen))
                        # - Simulate with crnsimulator python package (slower)
                        _odename = name + str(tlen)
                        _t0 = args.t0  if args.t0 > 0 and t_log else 1e-6
                        tfile = DiGraphSimulator(CG, _fname, nlist, p0, _t0, _t8,
                                                 t_lin=t_lin,
                                                 t_log=t_log,
                                                 jacobian=False,  # faster!
                                                 verb=(args.verbose > 0))
                        fail += 1

            except ExecError as e:
                # NOTE: This is a hack to avoid treekin simulations in the first place
                _odename = name + str(tlen)
                _t0 = args.t0  if args.t0 > 0 and t_log else 1e-6
                tfile = DiGraphSimulator(CG, _fname, nlist, p0, _t0, _t8,
                                         t_lin=t_lin,
                                         t_log=t_log,
                                         jacobian=False,  # faster!
                                         verb=(args.verbose > 1))

            # Get Results
            time_inc, iterations = CG.update_occupancies_tkn(tfile, nlist)
            # print time_inc, iterations


            if args.pyplot or dfh or args.stdout == 'drf':
                for rdata in CG.sorted_trajectories_iter(nlist, tfile, softmap):
                    [id_, tt_, oc_, ss_, en_] = rdata
                    if args.pyplot:
                        all_courses[id_].append((tt_, oc_))

                    fdata = "{} {} {} {:s} {:6.2f}\n".format(*rdata)
                    write_output(fdata, stdout=(args.stdout == 'drf'), fh = dfh)

            CG._total_time += time_inc

            # Prune
            dn, sr, rj = CG.prune(mocca=args.mocca, detailed=args.detailed_pruning)

            if args.draw_graphs:
                CG.to_json(_fname)

        if args.verbose:
            print("# Transcripton length: " + \
                "{}. Initial graph size: {}. Hidden graph size: {}".format(
                        tlen, len(nlist), len(CG)))
            print("#  Deleted {} nodes, {} still reachable, {} rejected deletions.".format(
                    dn, sr, rj))
            dtime = datetime.now()
            print("#  Computation time at current nucleotide: {} s".format(
                (dtime - atime).total_seconds()))

            fp_tot = ril.trafo.PROFILE['findpath-calls']
            fp_exp = ril.trafo.PROFILE['mfe'] + \
                     ril.trafo.PROFILE['hb'] + ril.trafo.PROFILE['feature']
            fp_cgr = ril.trafo.PROFILE['cogr']
            fp_prn = ril.trafo.PROFILE['prune']
            print("{} {} {} {} {} {} {}".format(tlen, len(nlist), 
                (dtime - atime).total_seconds(), fp_exp, fp_cgr, fp_prn, fp_tot))

            atime = dtime
            ril.trafo.PROFILE['findpath-calls'] = 0
            ril.trafo.PROFILE['mfe'] = 0
            ril.trafo.PROFILE['hb'] = 0
            ril.trafo.PROFILE['feature'] = 0
            ril.trafo.PROFILE['cogr'] = 0
            ril.trafo.PROFILE['prune'] = 0

    if args.verbose >= 1:
        print("# Treekin stats: {} default success, {} expo success, {} plusI success, {} fail".format( norm, expo, plusI, fail))

    # Write the last results
    if args.stdout == 'log' or lfh:
        fdata  = "Distribution of structures at the end:\n"
        fdata += "          {}\n".format(CG.transcript)
        for e, (ni, data) in enumerate(CG.sorted_nodes(), 1):
            sm = '-> {}'.format(CG.node[softmap[ni]]['identity']) if ni in softmap else ''
            fdata += "{:4d} {:4d} {} {:6.2f} {:6.4f} (ID = {:d} {:s})\n".format(
                tlen, e, ni[:tlen], data['energy'], data['occupancy'], data['identity'], sm)
        write_output(fdata, stdout=(args.stdout == 'log'), fh = lfh)

    # CLEANUP file handles
    if lfh: lfh.close()
    if dfh: dfh.close()

    # CLEANUP the /tmp/directory
    if not args.tmpdir:
        shutil.rmtree(_tmpdir)

    # Plot results
    if 'gr' in args.pyplot:
        plot_xmgrace(all_courses, filepath + '.gr')
        args.pyplot = filter(lambda x: x!='gr', args.pyplot)

    if args.pyplot:
        plot_simulation(all_courses, 
                steps = args.stop - args.start,
                t_ext = args.t_ext,
                t_end = args.t_end,
                fpath = filepath,
                formats = args.pyplot,
                title = args.name)

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        # formatter_class=argparse.RawTextHelpFormatter,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # formatter_class=argparse.MetavarTypeHelpFormatter,
        description='echo sequence | %(prog)s [options]')

    add_drtrafo_args(parser)

    args = parser.parse_args()

    main(args)
