#!/usr/bin/env python
#
# DrTransformer.py -- cotranscriptional folding.
#
import logging
import os
import sys
import math
import shutil
import tempfile
import argparse
from packaging import version

import RNA
import ribolands as ril
from ribolands.utils import parse_vienna_stdin
from ribolands.syswraps import sys_treekin
from ribolands.syswraps import SubprocessError, ExecError, check_version, VersionError
from ribolands.trafo import TrafoLandscape

def restricted_float(x):
    y = float(x)
    if y < 0.0 or y > 1.0:
        raise argparse.ArgumentTypeError(f"{x} not in range (0.0, 1.0)")
    return y

def sorted_trajectories(nlist, tfile, softmap = None):
    """ Yields the time course information using a treekin output file.

    Args:
      nlist (list): a list of nodes sorted by their energy
      tfile (str): treekin-output file name.
      softmap (dict, optional): A mapping to transfer occupancy between
        states. Likely not the most efficient implementation.

    Yields:
      list: ID, time, occupancy, structure, energy
    """
    with open(tfile) as tkn:
        for line in tkn:
            if line[0] == '#':
                continue
            course = list(map(float, line.strip().split()))
            time = course[0]
            if softmap:
                # Combine all results that belong together.
                raise NotImplementedError('No support for softmap atm.')
            for e, occu in enumerate(course[1:]):
                # is it above visibility threshold?
                yield time, nlist[e], occu 
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

def parse_drtrafo_args(parser):
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
    output.add_argument("--name", default = '', metavar = '<str>',
            help = """Name your output files, name the header of your plots, etc.
            this option overwrites the fasta-header.""")

    output.add_argument("--stdout", default = None, action = 'store',
            choices = ('log', 'drf', 'OFF'),
            help = """Choose STDOUT formats to follow the cotranscriptional
            folding progress in real time: *log*: a human readable output
            format.  *drf*: DrForna visalization input format. *OFF*: actively
            suppress output. The default (None) switches between *OFF*, if --logfile is
            specified, or *log* otherwise.""")

    output.add_argument("--logfile", action = "store_true",
            help = """Write verbose information to a file:
            {--outdir}/{--name}.log""")

    output.add_argument("--tmpdir", default = '', action = 'store', metavar = '<str>',
            help = """Specify path for storing temporary output files. These
            files will not be removed when the program terminates. Defaults to
            the default directory for temporary files on your System.""")

    output.add_argument("--outdir", default = '', action = 'store', metavar = '<str>',
            help = """Place regular output files, into this directory. Creates
            the directory if it does not exist. """)

    output.add_argument("--no-timecourse", action = "store_true",
            help = """Do not produce the time-course file (outdir/name.drf).""")

    output.add_argument("--t-inc", type = float, default = 1.2, metavar = '<flt>',
            help = """Adjust the plotting time resolution via the time-increment of
            the solver (t1 * t-inc = t2).""")

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
    algo.add_argument("--o-min", type = restricted_float, default = 0.1, metavar = '<flt>',
            help = """Occupancy threshold to determine which structures are
            relevant when transcribing a new nucleotide. A structure with
            occupancy o <  o-min / n gets pruned from the energy landscape, 
            where n is the number of simulated structures.
            """)

    algo.add_argument("--t-fast", type = float, default = 0.001, metavar = '<flt>',
            help = """Folding times faster than --t-fast are considered
            instantaneous.  Structural transitions that are faster than
            --t-fast are considerd part of the same macrostate. Directly
            translates to an energy barrier separating conforations using:
            dG = -RT*ln((1/t-fast)/k0). None: t-fast = 1/k_0 """)

    algo.add_argument("--minh", type = float, default = None, metavar = '<flt>',
            # An alternative to specify --t-fast in terms of a barrier height.
            help = argparse.SUPPRESS)

    algo.add_argument("--t-slow", type = float, default = 360000, metavar = '<flt>',
            help = """Only accept new structures as neighboring conformations if
            the transition is faster than --t-slow. This parameter may be
            useful to prevent numeric instabilities, otherwise, better avoid
            it.""")

    algo.add_argument("--maxh", type = float, default = None, metavar = '<flt>',
            # An alternative to specify --t-fast in terms of a barrier height.
            help = argparse.SUPPRESS)

    algo.add_argument("--force", action = "store_true", 
            # Enforce a setting against all warnings.
            help = argparse.SUPPRESS)

    # NOTE: findpath width is flexible by default, easy to implement though.
    #algo.add_argument("--findpath-search-width", type = int, default = 0, metavar = '<int>',
    #        help = """Search width for the *findpath* heuristic. Higher values
    #        increase the chances to find energetically better transition state
    #        energies.""")

    # NOTE: there is only one mode available right now.
    #algo.add_argument('--structure-search-mode', default = 'default',
    #        choices = ('default', 'fullconnect'),
    #        help = argparse.SUPPRESS)
    #        #help = """Specify one of three modes: *default*: find new secondary
    #        #structures using both the current MFE structure and fraying
    #        #neighbors.  *mfe-only*: only find the current MFE structure at
    #        #every transcription step.  *fraying-only*: only find local
    #        #fraying neighbors at every transcription step.""")

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
    return

def write_output(data, stdout = False, fh = None):
    # Helper function to print data to filehandle, STDOUT, or both.
    if stdout:
        sys.stdout.write(data)
    if fh:
        fh.write(data)
    return

def set_handle_verbosity(h, v):
    if v == 0:
        h.setLevel(logging.WARNING)
    elif v == 1:
        h.setLevel(logging.INFO)
    elif v == 2:
        h.setLevel(logging.DEBUG)
    elif v >= 3:
        h.setLevel(logging.NOTSET)

def main():
    """ DrTransformer - cotranscriptional folding. 
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'DrTransformer: RNA folding kinetics during transcription.')
    parse_drtrafo_args(parser)
    parse_model_details(parser)
    args = parser.parse_args()

    # ~~~~~~~~~~~~~
    # Logging Setup 
    # ~~~~~~~~~~~~~
    title = 'DrTransformer: RNA folding kinetics during transcription.'
    logger = logging.getLogger('ribolands')
    logger.setLevel(logging.DEBUG)

    banner = "{} {}".format(title, ril.__version__)
    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    set_handle_verbosity(ch, args.verbose)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info(banner)

    (name, fullseq) = parse_vienna_stdin(sys.stdin, chars='ACGUNacgun')

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
    else:
        filepath = args.name

    dfh = open(filepath + '.drf', 'w') if not args.no_timecourse else None
    lfh = open(filepath + '.log', 'w') if args.logfile else None

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

    if args.minh:
        logger.warning('Overwriting t-fast parameter.')
        args.t_fast = 1/(args.k0 * math.exp(-args.minh/_RT))
    else:
        args.minh = max(0, -_RT * math.log(1 / args.t_fast / args.k0))
    logger.info(f'--t-fast: {args.t_fast} s => {args.minh} kcal/mol barrier height ' + 
                f'and {1/args.t_fast} /s rate at k0 = {args.k0}')

    if args.maxh:
        logger.warning('Overwriting t-slow parameter.')
        args.t_slow = 1/(args.k0 * math.exp(-args.maxh/_RT))
    else:
        args.maxh = -_RT * math.log(1 / args.t_slow / args.k0)
    logger.info(f'--t-slow {args.t_slow} s => {args.maxh} kcal/mol barrier height ' +
                f'and {1/args.t_slow} /s rate at k0 = {args.k0}')

    if not args.force and args.t_fast and args.t_fast * 10 > args.t_ext:
        raise SystemExit('ERROR: Conflicting Settings: ' + 
                'Arguments must be such that "--t-fast" * 10 > "--t-ext".\n' + 
                '       => An instant folding time must be at least 10x shorter than ' +
                'the time of nucleotide extension. You may use --force to ignore this setting.')

    if not args.force and args.t_slow and args.t_end * 100 > args.t_slow:
        raise SystemExit('ERROR: Conflicting Settings: ' + 
                'Arguments must be such that "--t-slow" < 100 * "--t-end".\n' + 
                '       => A negligibly long folding time must be at least 100x longer than ' +
                'the final simulation time. You may use --force to ignore this error.')

    try:
        check_version(args.treekin, ril._MIN_TREEKIN_VERSION)
    except ExecError:
        pass

    if version.parse(RNA.__version__) < version.parse(ril._MIN_VIENNARNA_VERSION):
        raise VersionError('ViennaRNA', RNA.__version__, ril._MIN_VIENNARNA_VERSION)

    ############################
    # ~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Start with DrTransformer #
    # ~~~~~~~~~~~~~~~~~~~~~~~~ #
    ############################

    if args.paramFile:
        RNA.read_parameter_file(args.paramFile)

    # Set model details.
    vrna_md = RNA.md()
    vrna_md.noLP = 1
    vrna_md.logML = 0
    vrna_md.temperature = args.temp
    vrna_md.dangles = args.dangles
    vrna_md.special_hp = not args.noTetra
    vrna_md.noGU = args.noGU
    vrna_md.noGUclosure = args.noClosingGU
    #vrna_md.gquad = 0 G-Quads cannot be turned on.


    # Write logging output
    if args.stdout == 'log' or lfh:
        fdata  = "# File generated using DrTransformer v{}\n".format(ril.__version__)
        fdata += "#\n"
        fdata += "# >{}\n# {} \n".format(name, fullseq)
        fdata += "#\n"
        fdata += "# Co-transcriptional folding parameters:\n"
        fdata += "# --t-ext: {} sec\n".format(args.t_ext)
        fdata += "# --t-end: {} sec\n".format(args.t_end)
        fdata += "# --start: {}\n".format(args.start)
        fdata += "# --stop: {}\n".format(args.stop)
        fdata += "#\n"
        fdata += "# Algorithm parameters:\n"
        fdata += "# --o-min: {}\n".format(args.o_min)
        fdata += "# --t-fast: {} sec\n".format(args.t_fast)
        fdata += "# --t-slow: {} sec\n".format(args.t_slow)
        #fdata += "# --findpath-search-width: {}\n".format(args.findpath_search_width)
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
        fdata += "# Tanscription Step | Energy-sorted structure index | Structure | Energy "
        fdata += "| [Occupancy-t0 Occupancy-t8] | Structure ID (-> Plotting ID)\n"
        write_output(fdata, stdout=(args.stdout == 'log'), fh = lfh)

    # Write DrForna output
    if args.stdout == 'drf' or dfh:
        # Dictionary to store time course data.
        all_courses = dict()
        fdata = "id time conc struct energy\n"
        write_output(fdata, stdout = (args.stdout == 'drf'), fh = dfh)

    # initialize a directed conformation graph
    TL = TrafoLandscape(fullseq, vrna_md)
    TL.k0 = args.k0
    TL.minh = int(round(args.minh*100)) if args.minh else 0
    TL.maxh = int(round(args.maxh*100)) if args.maxh else 0
    TL._transcript_length = args.start - 1

    tprofile = []
    psites = dict()
    if args.pause_sites:
        for term in args.pause_sites:
            site, pause = term.split('=')
            psites[int(site)] = float(pause)

    #############
    # ~~~~~~~~~ #
    # Main loop #
    # ~~~~~~~~~ #
    #############
    for tlen in range(args.start, args.stop):
        time = TL.total_time
        # Get new nodes and connect them.
        logger.info(f'Expanding from {len(list(TL.active_local_mins))} local minima. ' + 
                    f'{len(list(TL.active_nodes))=} {len(list(TL.inactive_nodes))=}.')

        nn = TL.expand(mfree = args.min_fraying)
        logger.info(f'Found {len(nn)} new nodes during expansion.')
        cn, ce = TL.get_coarse_network()
        logger.info(f'Simulation network size: nodes = {cn}, edges = {ce}.')

        # ~~~~~~~~ #
        # Simulate #
        # ~~~~~~~~ #
        _fname = _tmpdir + '/' + name + '-' + str(tlen)
        _t0 = args.t0
        if tlen == args.stop - 1:
            _t8 = args.t_end 
        else: 
            _t8 = psites[tlen] if tlen in psites else args.t_ext
        tprofile.append(_t8)
        # Fallback mode: Python solver.
        (t_lin, t_log) = (None, args.t_log) if tlen == args.stop - 1 else (args.t_lin, None)

        # Produce input for treekin simulation
        nlist, bbfile, brfile, bofile, p0 = TL.get_simulation_files_tkn(_fname) 
        # Store the occupancy before the simulation for logfile output.
        odict = {n: TL.nodes[n]['occupancy'] for n in nlist}

        if len(nlist) == 1: # Fake simulation.
            node = nlist[0]
            ne = TL.nodes[node]['energy']/100
            no = TL.nodes[node]['occupancy']
            ni = TL.nodes[node]['identity']
            ss = node[:tlen]
            if args.stdout == 'log' or lfh:
                fdata = f"{tlen:4d} {1:4d} {ss} {ne:6.2f} [{no:6.4f} -> {no:6.4f}] ID = {ni}\n"
                write_output(fdata, stdout = (args.stdout == 'log'), fh = lfh)
            if args.stdout == 'drf' or dfh:
                #if ni not in all_courses:
                #    all_courses[ni] = [(time, 0)]
                #    fdata = f"{ni} {time:03.9f} {0:03.4f} {ss} {ne:6.2f}\n"
                #    write_output(fdata, stdout = (args.stdout == 'drf'), fh = dfh)
                all_courses[ni] = [(time + _t8, 1)]
                fdata = f"{ni} {time+_t8:03.9f} {no:03.4f} {ss} {ne:6.2f}\n"
                write_output(fdata, stdout = (args.stdout == 'drf'), fh = dfh)
            TL.total_time += _t8
            continue
        try:  # - Simulate with treekin
            tfile, _ = sys_treekin(_fname, bbfile,
                    treekin = args.treekin, p0 = p0, t0 = _t0, ti = args.t_inc, t8 = _t8,
                    binrates = True, mpack_method = args.mpack_method, quiet = True,
                    exponent = False, useplusI = False, force = True)
        except SubprocessError:
            try:  # - Simulate with treekin and --exponent
                tfile, _ = sys_treekin(_fname, bbfile,
                        treekin = args.treekin, p0 = p0, t0 = _t0, ti = args.t_inc, t8 = _t8,
                        binrates = True, mpack_method = args.mpack_method, quiet = True,
                        exponent = True, useplusI = False, force = True)
            except SubprocessError:
                try:  # - Simulate with treekin and --useplusI
                    tfile, _ = sys_treekin(_fname, bbfile,
                            treekin = args.treekin, p0 = p0, t0 = _t0, ti = args.t_inc, t8 = _t8,
                            binrates = True, mpack_method = args.mpack_method, quiet = True,
                            exponent = True, useplusI = True, force = True)
                except SubprocessError:
                    raise SystemExit(f"After {tlen} nucleotides: treekin cannot find a solution!")

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # Process simulation results #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # Update occupancies from treekin results. 
        time_inc = TL.update_occupancies_tkn(tfile, nlist)

        # Print the current state *after* the simulation.
        if args.stdout == 'log' or lfh:
            for e, node in enumerate(nlist, 1):
                ne = TL.nodes[node]['energy']/100
                po = odict[node]
                no = TL.nodes[node]['occupancy']
                ni = TL.nodes[node]['identity']
                fdata = f"{tlen:4d} {e:4d} {node[:tlen]} {ne:6.2f} [{po:6.4f} -> {no:6.4f}] ID = {ni}\n"
                write_output(fdata, stdout = (args.stdout == 'log'), fh = lfh)

        if args.stdout == 'drf' or dfh:
            lt = time
            for (stime, node, occu) in sorted_trajectories(nlist, tfile):
                tt = time + stime
                ni = TL.nodes[node]['identity']
                if occu < 0.001: 
                    occu = 0
                if occu == 0 and ni not in all_courses:
                    continue
                ne = TL.nodes[node]['energy']/100
                #if ni not in all_courses:
                #    if lt < tt:
                #        all_courses[ni] = [(lt, 0)]
                #        fdata = f"{ni} {lt:03.9f} {0:03.4f} {node[:tlen]} {ne:6.2f}\n"
                #        write_output(fdata, stdout = (args.stdout == 'drf'), fh = dfh)
                #    else:
                #        all_courses[ni] = []
                all_courses[ni] = [(tt, occu)]
                fdata = f"{ni} {tt:03.9f} {occu:03.4f} {node[:tlen]} {ne:6.2f}\n"
                write_output(fdata, stdout = (args.stdout == 'drf'), fh = dfh)
                lt = tt

        # ~~~~~ #
        # Prune #
        # ~~~~~ #
        if args.o_min > 0:
            # Adjust p-min to size of current structure space.
            pmin = args.o_min / len(nlist)
            logger.info(f'Pruning nodes with occupancy below {pmin}.')
            TL.prune(pmin)
        if args.draw_graphs:
            raise NotImplementedError('Cannot draw graphs right now.')
            TL.graph_to_json(_fname)
        TL.total_time += time_inc

    # Write the last results
    if args.stdout == 'log' or lfh:
        fdata  = "# Distribution of structures at the end:\n"
        fdata += "#         {}\n".format(TL.transcript)
        for e, node in enumerate(sorted(TL.active_local_mins, key = lambda x: TL.nodes[x]['energy']), 1):
            ne = TL.nodes[node]['energy']/100
            no = TL.nodes[node]['occupancy']
            ni = TL.nodes[node]['identity']
            fdata += f"{tlen:4d} {e:4d} {node[:tlen]} {ne:6.2f} {no:6.4f} ID = {ni}\n"
        write_output(fdata, stdout=(args.stdout == 'log'), fh = lfh)
    logger.info(f"Transkription profile: {tprofile}\n")

    # CLEANUP file handles
    if lfh: lfh.close()
    if dfh: dfh.close()
    # CLEANUP the /tmp/directory
    if not args.tmpdir:
        shutil.rmtree(_tmpdir)
    return

if __name__ == '__main__':
    main()
