#!/usr/bin/env python

#
# DrTransformer.py -- cotranscriptional folding.
#
# written by Stefan Badelt (stef@tbi.univie.ac.at)
#

from __future__ import division, print_function
from builtins import map
from builtins import zip
from builtins import str
from builtins import range

import os
import sys
import math
import shutil
import tempfile
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from datetime import datetime
import seaborn

import RNA
import ribolands as ril
from ribolands.syswraps import SubprocessError, ExecError, check_version, VersionError
from ribolands.crnwrapper import DiGraphSimulator
from ribolands.trafo import TrafoLandscape

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

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
        tprofile,
        fpath,
        formats = ['pdf'],
        title = ''):
    """
    """
    seaborn.set_style("darkgrid")

    t_end = tprofile[-1]

    # Get the relevant arguments from args
    lin_time = sum(tprofile[:-1])
    log_time = lin_time + t_end if (lin_time + t_end) >= lin_time * 10 else lin_time * 10

    # Do the plotting
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['right'].set_visible(False)
    ax.set_ylim([-0.05, 1.05])
    ax.set_xscale('linear')

    # Make the second part of the plot logarithmic
    offset = 0.00001
    ax.set_xlim((0, lin_time + offset))
    divider = make_axes_locatable(ax)
    axLog = divider.append_axes("right", size=2.5, pad=0, sharey=ax)
    axLog.set_xscale('log')
    axLog.set_xlim((lin_time + offset, log_time))
    axLog.set_ylim([-0.05, 1.05])
    axLog.yaxis.set_visible(False)
    axLog.spines['left'].set_visible(False)

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

    #for tlen in range(len(tprofile)-1):
    #    ax.axvline(x = sum(tprofile[:tlen]), linewidth = 0.01, color = 'black', linestyle = '--')

    # Add ticks for 1 minute, 1 hour, 1 day, 1 year
    axLog.axvline(x = lin_time, linewidth = 3, color = 'black', linestyle = '-') 
    axLog.axvline(x = 60, linewidth = 1, color = 'black', linestyle = '--') # 1 minute
    axLog.axvline(x = 3600, linewidth = 1, color = 'black', linestyle = '--') # 1 hour
    axLog.axvline(x = 86400, linewidth = 1, color = 'black', linestyle = '--') # 1 day
    axLog.axvline(x = 31536000, linewidth = 1, color = 'black', linestyle = '--') # 1 year
    plt.legend()

    ax.set_ylabel('occupancy', fontsize = 11)
    ax.set_xlabel('time [s]', ha = 'center', va = 'center', fontsize = 11)
    ax.xaxis.set_label_coords(.9, -0.15)

    for ending in formats:
        pfile = fpath + '.' + ending
        plt.savefig(pfile, bbox_inches = 'tight')

    return

def parse_model_details(parser):
    """ ViennaRNA Model Details Argument Parser.  """
    model = parser.add_argument_group('ViennaRNA model details')

    model.add_argument("-T", "--temp", type = float, default = 37.0, metavar = '<flt>',
        help = 'Rescale energy parameters to a temperature of temp C.')

    model.add_argument("-4", "--noTetra", action = "store_true",
        help = 'Do not include special tabulated stabilizing energies for tri-, tetra- and hexaloop hairpins.')

    model.add_argument("-d", "--dangles", type = int, default = 2, metavar = '<int>',
        help = 'How to treat "dangling end" energies for bases adjacent to helices in free ends and multi-loops.')

    model.add_argument("--noGU", action = "store_true",
        help = 'Do not allow GU/GT pairs.')

    model.add_argument("--noClosingGU", action = "store_true",
        help = 'Do not allow GU/GT pairs at the end of helices.')

    model.add_argument("-P", "--paramFile", action = "store", default = None, metavar = '<str>',
        help = 'Read energy parameters from paramfile, instead of using the default parameter set.')

def add_drtrafo_args(parser):
    """ A collection of arguments that are used by DrTransformer """

    environ = parser.add_argument_group('DrTransformer dependencies')
    output = parser.add_argument_group('DrTransformer output')
    trans = parser.add_argument_group('Transcription parameters')
    algo  = parser.add_argument_group('DrTransformer algorithm')

    ###################
    # Default options #
    ###################
    parser.add_argument('--version', action = 'version', 
            version = '%(prog)s ' + ril.__version__)
    parser.add_argument("-v", "--verbose", action = 'count', default = 0,
            help = """Track process by writing verbose output to STDOUT during
            calculations. Use --logfile if you want to see *just* verbose
            information via STDOUT.  The verbose output can be visualized using
            the script DrProfile.py. """)

    ##############################
    # DrTransformer dependencies #
    ##############################
    environ.add_argument("--treekin", default = 'treekin', action = 'store', metavar = '<str>', 
            help = "Path to the *treekin* executable.")

    environ.add_argument("--mpack-method", action = 'store',
            choices = (['FLOAT128', 'LD', 'QD', 'DD']),
            help = """Increase the precision of treekin simulations. Requires a development
                    version of treekin with mpack support.""")

    ########################
    # DrTransformer output #
    ########################
    output.add_argument("--stdout", default = None, action = 'store',
            choices = ('log', 'drf', 'OFF'),
            help = """Choose STDOUT formats to follow the cotranscriptional
            folding progress in real time: *log*: a human readable output
            format.  *drf*: DrForna visalization input format. *OFF*: actively
            suppress output. The default (None) switches between *OFF*, if --logfile is
            specified, or *log* otherwise.""")

    output.add_argument("--tmpdir", default = '', action = 'store', metavar = '<str>',
            help = """Specify path for storing temporary output files. These
            files will not be removed when the program terminates. Defaults to
            the default directory for temporary files on your System.""")

    output.add_argument("--outdir", default = '', action = 'store', metavar = '<str>',
            help = """Place regular output files, into this directory. Creates
            the directory if it does not exist. """)

    output.add_argument("--name", default = '', metavar = '<str>',
            help = """Name your output files, name the header of your plots, etc.
            this option overwrites the fasta-header.""")

    output.add_argument("--logfile", action = "store_true",
            help = """Write verbose information to a file:
            {--outdir}/{--name}.log""")

    output.add_argument("--visualize", nargs = '+', default = '', action = 'store',
            choices = ('pdf', 'svg', 'png', 'gr', 'drf'),
            help = """Plot the simulation using matplotlib (pdf, svg, png) and/or
            write an input file for xmgrace (gr) and/or write an input file for
            DrForna (drf). Interpret the legend using STDOUT or --logfile.
            Files: {--outdir}/{--name}.{--visualize}""")

    output.add_argument("--t-inc", type = float, default = 1.2, metavar = '<flt>',
            help = """Adjust the plotting time resolution via the time-increment of
            the solver (t1 * t-inc = t2).""")

    output.add_argument("--plot-minh", type = float, default = 0, metavar = '<flt>',
            help = """Reduce the resolution of visualized structures. This
            parameter merges structures which are separated by a barrier
            smaller than --plot-minh *for visualzation only*. The dynamics will
            still be caculated based on the more detailed network. This
            parameter can only be effectve if it leads to a stronger
            coarse-graining than the simulation paramter --t-fast.  
            
            This parameter corresponds to folding time as: 
            t = 1/(k0*exp(-dG/RT)), where dG is the parameter --plot-minh.
            E.g.  --plot-minh: 5.1 kcal/mol ~=  0.02 s [default t-ext]. """)

    output.add_argument("--draw-graphs", action = "store_true",
            #help = """Export every landscape as json file. Uses --tempdir. """)
            help = argparse.SUPPRESS)

    output.add_argument("--t-lin", type = int, default = 30, metavar = '<int>',
            #help = """Evenly space output *t-lin* times during transcription on a linear time scale.""")
            help = argparse.SUPPRESS)
    output.add_argument("--t-log", type = int, default = 300, metavar = '<int>',
            #help = """Evenly space output *t-log* times after transcription on a logarithmic time scale.""")
            help = argparse.SUPPRESS)

    output.add_argument("--t0", type = float, default = 0, metavar = '<flt>',
            help = argparse.SUPPRESS)

    ############################
    # Transcription parameters #
    ############################
    trans.add_argument("--t-ext", type = float, default = 0.02, metavar = '<flt>',
            help = """Transcription speed, i.e. time per nucleotide extension
            [seconds per nucleotide].""")
    trans.add_argument("--t-end", type = float, default = 60, metavar = '<flt>',
            help = "Post-transcriptional simulation time [seconds].")
    trans.add_argument("--pause-sites", nargs='+', metavar='<int>=<flt>',
            help="""Transcriptional pausing sites.  E.g. \"--pause-sites 82=2e3
            84=33\" alters the simulation time at nucleotides 82 and 84 to 2000
            and 33 seconds, respectively. """)
    trans.add_argument("--start", type = int, default = 1, metavar = '<int>',
            help = "Start transcription at this nucleotide.")
    trans.add_argument("--stop", type = int, default = None, metavar = '<int>',
            help = "Stop transcription before this nucleotide")

    ###########################
    # DrTransformer algorithm #
    ###########################
    algo.add_argument("--track-basins", type = float, default = None, metavar = '<flt>',
            # Specify a number of basins to follow, and adjust coares-graining accordingly.
            help = argparse.SUPPRESS)

    algo.add_argument("--o-min", type = restricted_float, default = 0.1, metavar = '<flt>',
            help = """Occupancy threshold to determine which structures are
            relevant when transcribing a new nucleotide. A structure with
            occupancy o <  o-min / n gets pruned from the energy landscape, 
            where n is the number of simulated structures.
            """)

    algo.add_argument("--t-fast", type = float, default = None, metavar = '<flt>',
            help = """Folding times faster than --t-fast are considered
            instantaneous.  Structural transitions that are faster than
            --t-fast are considerd part of the same macrostate. Directly
            translates to an energy barrier separating conforations using:
            dG = -RT*ln((1/t-fast)/k0). None: t-fast = 1/k_0 """)

    algo.add_argument("--minh", type = float, default = None, metavar = '<flt>',
            # An alternative to specify --t-fast in terms of a barrier height.
            help = argparse.SUPPRESS)

    algo.add_argument("--force", action = "store_true", 
            # Enforce a setting against all warnings.
            help = argparse.SUPPRESS)

    algo.add_argument("--t-slow", type = float, default = None, metavar = '<flt>',
            help = """Only accept new structures as neighboring conformations if
            the transition is faster than --t-slow. This parameter may be
            useful to prevent numeric instabilities, otherwise, better avoid
            it.""")

    algo.add_argument("--findpath-search-width", type = int, default = 20, metavar = '<int>',
            help = """Search width for the *findpath* heuristic. Higher values
            increase the chances to find energetically better transition state
            energies.""")

    algo.add_argument('--structure-search-mode', default = 'default',
            choices = ('default', 'fullconnect'),
            help = argparse.SUPPRESS)
            #help = """Specify one of three modes: *default*: find new secondary
            #structures using both the current MFE structure and fraying
            #neighbors.  *mfe-only*: only find the current MFE structure at
            #every transcription step.  *fraying-only*: only find local
            #fraying neighbors at every transcription step.""")

    algo.add_argument("--min-fraying", type = int, default = 6, metavar = '<int>',
            help = """Minimum number of freed bases during helix fraying.
            Fraying helices can vary greatly in length, starting with at
            least two base-pairs. This parameter defines the minimum amount of
            bases freed by helix fraying. For example, 6 corresponds to a
            stack of two base-pairs and a loop region of 2 nucleotides. If less
            bases are freed and there exists a nested stacked helix, this helix
            is considered to fray as well.""")

    algo.add_argument("--k0", type = float, default = 2e5, metavar = '<flt>',
            help = """Arrhenius rate constant (pre-exponential factor). Adjust
            this constant of the Arrhenius equation to relate free energy
            changes to experimentally determined folding time [seconds].""")

    # DEPRECATED
    algo.add_argument("--min-occupancy", type = restricted_float, default = None, metavar = '<flt>',
            help = argparse.SUPPRESS)

    output.add_argument("--soft-minh", type = float, default = 0, metavar = '<flt>',
            #NOTE: Deprecated: use --plot-minh!
            help = argparse.SUPPRESS)

    output.add_argument("--tinc", type = float, default = None, metavar = '<flt>',
            #NOTE: Deprecated: use --t-inc!
            help = argparse.SUPPRESS)

    # TESTING
    algo.add_argument("--fake-pruning", type = float, default = None, metavar = '<flt>',
            # An alternative to specify --t-fast in terms of a barrier height.
            help = argparse.SUPPRESS)

    return


def write_output(data, stdout = False, fh = None):
    # Helper function to print data to filehandle, STDOUT, or both.
    if stdout:
        sys.stdout.write(data)
    if fh:
        fh.write(data)
    return

def main(args):
    """ DrTransformer - cotranscriptional folding """
    (name, fullseq) = ril.parse_vienna_stdin(sys.stdin, chars='ACGUNacgun')

    # Argument deprecation Warnings
    if args.tinc:
        print('DEPRECATION WARNING: Argument --tinc is deprecated, please use --t-inc.')
        args.t_inc = args.tinc

    if args.soft_minh:
        print('DEPRECATION WARNING: Argument --soft-minh is deprecated, please use --plot-minh.')
        args.plot_minh = args.soft_minh

    if args.min_occupancy:
        print('DEPRECATION WARNING: Argument --min-occupancy is deprecated, please use --o-min.')
        args.o_min = args.min_occupancy

    # Adjust arguments, prepare simulation
    if args.name == '':
        args.name = name
    else:
        name = args.name

    if os.path.split(args.name)[0]:
        raise SystemExit('ERROR: Argument "--name" must not contain filepath.')

    if args.outdir:
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        filepath = args.outdir + '/' + args.name
    else :
        filepath = args.name

    dfh = open(filepath + '.drf', 'w') if 'drf' in args.visualize else None
    lfh = open(filepath + '.log', 'w') if args.logfile else None

    args.pyplot = [x for x in args.visualize if x not in ('drf')]

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
    if args.temp != 37.0:
        kelvin = 273.15 + args.temp
        _RT = (_RT / 310.15) * kelvin

    if args.stop is None:
        args.stop = len(fullseq) + 1

    if args.t_end < args.t_ext:
        raise SystemExit('ERROR: Conflicting Settings: ' + \
                'Arguments must be such that "--t-end" >= "--t-ext"')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Adjust the simulation window for treekin:
    #
    #   k0
    #   2e5  atu/s    50 nuc/s               1e-inf /s
    #   4000 /ext     1  nuc/ext             1e-inf /ext
    #   |------$------|-------------------$----------->
    #          k-fast                     k-slow
    #          --minh                     --maxh
    #          |------|---------------|--->
    #          t0     t_ext 0.02 s    t_end = 86400 s
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
    # post-transcriptional simulation time --t-end
    #
    # Parameters:
    # k0 = maximum folding rate /s
    # t-ext = time for chain elongation
    # t-end = post-transcriptional simulation time
    # t0 = first simulation output time (<< t8)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    if args.t_fast is None:
        if args.minh:
            args.t_fast = 1 / (args.k0 * math.exp(-args.minh/_RT))
        else:
            args.t_fast = 1 / args.k0
    elif args.minh:
        raise SystemExit('ERROR: Conflicting settings: "--minh" and "--t-fast" must not be specified at the same time.')

    # Compare --t-fast, --t-slow, --t-ext and --t-end :
    if args.verbose:
        dG_min = -_RT * math.log(1 / args.t_fast / args.k0)
        print('# --t-fast: {} s => {} kcal/mol barrier height and {} /s rate at k0 = {}'.format(
            args.t_fast, dG_min, 1/args.t_fast, args.k0))
        if args.t_slow is not None:
            dG_max = -_RT * math.log(1 / args.t_slow / args.k0)
            print('# --t-slow {} s => {} kcal/mol barrier height and {} /s rate at k0 = {}'.format(
                args.t_slow, dG_max, 1/args.t_slow, args.k0))

    if not args.force and args.t_fast * 10 > args.t_ext:
        raise SystemExit('ERROR: Conflicting Settings: ' + \
                'Arguments must be such that "--t-fast" * 10 > "--t-ext"\n' + \
                'The rate for equilibration must be at least 10x faster than for nucleotide extension. You may use --force to ignore this error.')

    if not args.force and args.t_slow is not None and args.t_end * 100 > args.t_slow:
        raise SystemExit('ERROR: Conflicting Settings: ' + \
                'Arguments must be such that "--t-slow" < "--t-end" * 100\n' + \
                'A negligible rate must be at least 100x slower than the final simulation time --t-end. You may use --force to ignore this error.')

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

    if args.paramFile:
        RNA.read_parameter_file(args.paramFile)

    # Set model details.
    vrna_md = RNA.md()
    vrna_md.noLP = 1
    vrna_md.temperature = args.temp
    vrna_md.dangles = args.dangles
    vrna_md.logML = 0
    vrna_md.special_hp = not args.noTetra
    vrna_md.noGU = args.noGU
    vrna_md.noGUclosure = args.noClosingGU
    #vrna_md.gquad = 0 G-Quads cannot be turned on.

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
        fdata += "# --start: {}\n".format(args.start)
        fdata += "# --stop: {}\n".format(args.stop)
        fdata += "#\n"
        fdata += "# Algorithm paramters:\n"
        fdata += "# --o-min: {}\n".format(args.o_min)
        fdata += "# --t-fast: {} sec\n".format(args.t_fast)
        fdata += "# --t-slow: {} sec\n".format(args.t_slow)
        fdata += "# --findpath-search-width: {}\n".format(args.findpath_search_width)
        fdata += "# --min-fraying: {} nuc\n".format(args.min_fraying)
        fdata += "# --k0: {}\n".format(args.k0)
        fdata += "#\n"
        fdata += "# ViennaRNA model details:\n"
        fdata += "# --temp: {} C\n".format(args.temp)
        fdata += "# --dangles: {}\n".format(args.dangles)
        fdata += "# --paramFile: {}\n".format(args.paramFile)
        fdata += "# --noTetra: {}\n".format(args.noTetra)
        fdata += "# --noGU: {}\n".format(args.noGU)
        fdata += "# --noClosingGU: {}\n".format(args.noClosingGU)
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
    CG._k0 = args.k0
    CG.t_fast = args.t_fast
    dG_min = CG._dG_min # currently: CG._dG_min is set via t-fast
    CG.t_slow = args.t_slow

    CG.findpath_search_width = args.findpath_search_width

    CG._transcript_length = args.start-1

    #import statprof
    #statprof.start()

    tprofile = []
    psites = dict()
    if args.pause_sites:
        for term in args.pause_sites:
            site, pause = term.split('=')
            psites[int(site)]=float(pause)

    if args.verbose:
        atime = datetime.now()

    # now lets start...
    tnorm, tplusI, texpo, tfail, tfake = 0, 0, 0, 0, 0
    for tlen in range(args.start, args.stop):
        nn = CG.expand(exp_mode = args.structure_search_mode)

        mn = CG.coarse_grain()
        #if args.verbose:
        #    print("# Merged {} nodes after expanding {} new nodes.".format(len(mn), nn))

        if args.pyplot:
            ttt = CG.total_time
            if ttt == 0:
                all_courses.extend([[] for i in range(nn)])
            else:
                all_courses.extend([[(ttt, 0)] for i in range(nn)])

        ############
        # Simulate #
        ############
        _fname = _tmpdir + '/' + name + '-' + str(tlen)
        _t0 = args.t0  # if args.t0 > 0 else 1e-6

        if tlen == args.stop - 1 :
            _t8 = args.t_end 
        else: 
            _t8 = psites[tlen] if tlen in psites else args.t_ext
        tprofile.append(_t8)

        (t_lin, t_log) = (None, args.t_log) if tlen == args.stop - \
            1 else (args.t_lin, None)

        # produce input for treekin simulation
        [bfile, rfile, p0, nlist] = CG.get_simulation_files_tkn(_fname)

        softmap = dict()
        if args.plot_minh and args.plot_minh > CG._dG_min:
            copyCG = CG.graph_copy()
            softmap = copyCG.coarse_grain(dG_min=args.plot_minh)
            del copyCG

        # Print the current state *before* the simulation starts.
        if args.stdout == 'log' or lfh:
            for e, (ni, data) in enumerate(nlist, 1):
                sm = '-> {}'.format([CG.nodes[tr]['identity'] for tr in softmap[ni]]) if ni in softmap else ''
                fdata = "{:4d} {:4d} {} {:6.2f} {:6.4f} ID = {:d} {:s}\n".format(
                    tlen, e, ni[:tlen], data['energy'], data['occupancy'], data['identity'], sm)
                write_output(fdata, stdout=(args.stdout == 'log'), fh = lfh)

        if args.verbose:
            itime = datetime.now()

        dn, sr, rj = 0, 0, 0
        if len(nlist) == 1:
            # Fake simulation results for DrForna
            tfake += 1
            CG.total_time += _t8
            if args.stdout == 'drf' or dfh:
                ss = nlist[0][0]
                fdata = "{:d} {:03.9f} {:03.4f} {:s} {:6.2f}\n".format(CG.nodes[ss]['identity'],
                        CG.total_time, 1.0, ss[:len(CG.transcript)], CG.nodes[ss]['energy'])
                write_output(fdata, stdout=(args.stdout == 'drf'), fh = dfh)

            if args.pyplot:
                ss = nlist[0][0]
                ident = CG.nodes[ss]['identity']
                all_courses[ident].append((CG.total_time, 1.0))

        else:
            seq = ''
            # sometimes bfile causes a segfault, so let's leave it out.
            bfile = None
            try:  # - Simulate with treekin
                tfile, _ = ril.sys_treekin_051(_fname, rfile,
                        treekin = args.treekin, p0 = p0, t0 = _t0, ti = args.t_inc, t8 = _t8,
                        binrates = True, mpack_method = args.mpack_method, quiet = True,
                        exponent = False, useplusI = False, force = True, verbose = (args.verbose > 1))
                tnorm += 1
            except SubprocessError:
                try:  # - Simulate with treekin and --exponent
                    tfile, _ = ril.sys_treekin_051(_fname, rfile,
                            treekin = args.treekin, p0 = p0, t0 = _t0, ti = args.t_inc, t8 = _t8,
                            binrates = True, mpack_method = args.mpack_method, quiet = True,
                            exponent = True, useplusI = False, force = True, verbose = (args.verbose > 1))
                    texpo += 1
                except SubprocessError:
                    try:  # - Simulate with treekin and --useplusI
                        tfile, _ = ril.sys_treekin_051(_fname, rfile,
                                treekin = args.treekin, p0 = p0, t0 = _t0, ti = args.t_inc, t8 = _t8,
                                binrates = True, mpack_method = args.mpack_method, quiet = True,
                                exponent = True, useplusI = True, force = True, verbose = (args.verbose > 1))
                        tplusI += 1
                    except SubprocessError:
                        if args.verbose > 1:
                            print("After {} nucleotides: treekin cannot find a solution!".format(tlen))
                        # - Simulate with crnsimulator python package (slower)
                        _odename = name + str(tlen)
                        _t0 = args.t0  if args.t0 > 0 and t_log else 1e-6
                        tfile = DiGraphSimulator(CG, _fname, nlist, p0, _t0, _t8,
                                                 t_lin = t_lin,
                                                 t_log = t_log,
                                                 jacobian = False,  # faster!
                                                 verb = (args.verbose > 1))
                        tfail += 1

            ## NOTE: Enable this hack to avoid treekin simulations in the first place
            #except ExecError as e:
            #    _odename = name + str(tlen)
            #    _t0 = args.t0  if args.t0 > 0 and t_log else 1e-6
            #    tfile = DiGraphSimulator(CG, _fname, nlist, p0, _t0, _t8,
            #                             t_lin = t_lin,
            #                             t_log = t_log,
            #                             jacobian = False,  # faster!
            #                             verb = (args.verbose > 1))

            # Get Results
            time_inc, iterations = CG.update_occupancies_tkn(tfile, nlist)

            if args.pyplot or dfh or args.stdout == 'drf':
                # Avoid printing the initial time point of sequence elongation,
                # it seems confusing and redundant.
                fdata, itime = None, None
                for rdata in CG.sorted_trajectories_iter(nlist, tfile, softmap):
                    [id_, tt_, oc_, ss_, en_] = rdata

                    if fdata is None:
                        if itime is None or itime == tt_:
                            itime = tt_
                            continue

                    # id time conc struct energy
                    fdata = "{:d} {:03.9f} {:03.4f} {:s} {:6.2f}\n".format(*rdata)

                    # NOTE: It seems reasonable to activate the following line
                    # to reduce the plot size. However, it doesn't help much
                    # and it might introduce a visualization bug for some weird
                    # examples, so let's leave it for now:
                    #if oc_ < 0.001: continue

                    if args.pyplot:
                        all_courses[id_].append((tt_, oc_))
                    write_output(fdata, stdout = (args.stdout == 'drf'), fh = dfh)



            CG.total_time += time_inc

            # Prune
            if args.o_min > 0:
                # adjust o-min to size of current structure space:
                pmin = args.o_min / len(nlist)
                dn, sr = CG.prune(pmin, detailed = True, keep_reachables = False)
            else:
                dn = sr = 0

            if args.track_basins:
                # add or substract a 0.1 kcal/mol plot-minh for every structure
                # too much or too little.
                approach = (len(nlist)-dn - args.track_basins) / 10
                CG._dG_min = max(dG_min, CG._dG_min + approach)

            if args.draw_graphs:
                CG.graph_to_json(_fname)

        if args.verbose:
            nZedges = len([a for (a, b, d) in CG.edges(data = True) if d['saddle'] != float('inf') and CG.nodes[a]['active'] and CG.nodes[b]['active']])
            print("# Transcripton length: {}. Active graph size: {}. Non-zero transition edges: {}.  Hidden graph size: {}. Number of Edges: {}".format(
                tlen, len(nlist), nZedges, len(CG), CG.number_of_edges()))
            stime = datetime.now()
            algotime = (itime - atime).total_seconds()
            simutime = (stime - itime).total_seconds()
            tot_time = (stime - atime).total_seconds()
            print("# Computation time at current nucleotide: algo: {} simu: {} total: {}".format(algotime, simutime, tot_time))
            print("# Deleted {} nodes, {} still reachable.".format(dn, sr))
            print("# Treekin stats: {} default success, {} expo success, {} plusI success, {} fail, {} fake".format(tnorm, texpo, tplusI, tfail, tfake))
            fp_tot = ril.trafo.PROFILE['findpath-calls']
            fp_fr1 = ril.trafo.PROFILE['fraying1']
            fp_fr2 = ril.trafo.PROFILE['fraying2']
            fp_mfe = ril.trafo.PROFILE['mfe']
            fp_con = ril.trafo.PROFILE['connect']
            fp_cgr = ril.trafo.PROFILE['cogr']
            fp_prn = ril.trafo.PROFILE['prune']
            print("# Findpath stats: {} fraying, {} mfe connect, {} triangle connect, {} coarse-grain, {} prune, {} total.".format(fp_fr1+fp_fr2, fp_mfe, fp_con, fp_cgr, fp_prn, fp_tot))
            print("# dG-min: {} dG-max: {}".format(CG._dG_min, CG._dG_max))
            print("{}  {} {} {} {}  {} {} {}  {} {}  {} {} {} {} {}  {} {} {} {} {} {}  {} {}".format(tlen, 
                len(nlist), nZedges, len(CG), CG.number_of_edges(),
                algotime, simutime, tot_time,
                dn, sr,
                tnorm, texpo, tplusI, tfail, tfake,
                fp_fr1+fp_fr2, fp_mfe, fp_con, fp_cgr, fp_prn, fp_tot, 
                CG._dG_min, CG._dG_max))

            atime = stime
            ril.trafo.PROFILE['findpath-calls'] = 0
            ril.trafo.PROFILE['fraying1'] = 0
            ril.trafo.PROFILE['fraying2'] = 0
            ril.trafo.PROFILE['mfe'] = 0
            ril.trafo.PROFILE['connect'] = 0
            ril.trafo.PROFILE['cogr'] = 0
            ril.trafo.PROFILE['prune'] = 0
            sys.stdout.flush()

    # Write the last results
    if args.stdout == 'log' or lfh:
        fdata  = "# Distribution of structures at the end:\n"
        fdata += "#         {}\n".format(CG.transcript)
        for e, (ni, data) in enumerate(CG.sorted_nodes(), 1):
            sm = '-> {}'.format([CG.nodes[tr]['identity'] for tr in softmap[ni]]) if ni in softmap else ''
            fdata += "{:4d} {:4d} {} {:6.2f} {:6.4f} ID = {:d} {:s}\n".format(
                tlen, e, ni[:tlen], data['energy'], data['occupancy'], data['identity'], sm)
        write_output(fdata, stdout=(args.stdout == 'log'), fh = lfh)


    #statprof.stop()
    #statprof.display()

    # CLEANUP file handles
    if lfh: lfh.close()
    if dfh: dfh.close()

    # CLEANUP the /tmp/directory
    if not args.tmpdir:
        shutil.rmtree(_tmpdir)

    # Plot results
    if 'gr' in args.pyplot:
        plot_xmgrace(all_courses, filepath + '.gr')
        args.pyplot = list(filter(lambda x: x!='gr', args.pyplot))

    if args.pyplot:
        plot_simulation(all_courses, 
                tprofile = tprofile,
                fpath = filepath,
                formats = args.pyplot,
                title = args.name)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        # formatter_class = argparse.RawTextHelpFormatter,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        # formatter_class = argparse.MetavarTypeHelpFormatter,
        description = 'echo sequence | %(prog)s [options]')

    add_drtrafo_args(parser)
    parse_model_details(parser)

    args = parser.parse_args()

    main(args)
