#
# ribolands.utils
# 
# All sorts of useful stuff.
#

import logging
rlog = logging.getLogger(__name__)

import re
import sys
from struct import pack, unpack, calcsize

class RiboUtilsError(Exception):
    pass

class ProgressBar(object):
    def __init__(self, clockmax):
        if clockmax < 1:
            raise RiboUtilsError("Input clockmax {clockmax} should be a positive integer.")

        self.clock = 0
        self.clockmax = clockmax
        print('# |Progress' + ' ' * (clockmax - 8) + '|')
        sys.stdout.write('# |')

    def inc(self):
        self.clock += 1
        if self.clock % 10 == 0:
            sys.stdout.write(str(int(self.clock/10)))
        elif self.clock % 5 == 0:
            sys.stdout.write(",")
        else:
            sys.stdout.write(".")

        sys.stdout.flush()
        if self.clock == self.clockmax:
            sys.stdout.write("| 100%\n")

def argparse_add_arguments(parser,
                           RNAsubopt=False, barriers=False, treekin=False,
                           noLP=False, temperature=False, circ=False,
                           tmpdir=False, name=False, force=False, verbose=False,
                           start=False, stop=False, k0=False, tX=False, cutoff=False):
    """Commonly used argparse arguments when using the ribolands library.

    Args:
      parser (:obj:`argparse.ArgumentParser()`): Standard object provided by argparse.
      ``flags`` (bool): selectively choose which arguments to include.

    Returns:
      Modifies the first argument: parser.
    """

    # ViennaRNA suboptimal structures:
    if RNAsubopt:
        parser.add_argument("--RNAsubopt", default='RNAsubopt', action='store',
                            metavar='<str>',
                            help="Specify path to your *RNAsubopt* executable")
        # TODO: old version allows '-e' and sets default=0
        parser.add_argument("--s-ener", type=float, default=None, metavar='<flt>',
                            help="Set the energy range for suboptimal structure computation." +
                            " This will overwrite the options --s-maxe and --s-maxn.")
        parser.add_argument("--s-maxe", type=float, default=20, metavar='<flt>',
                            help="Set a the maximum subopt range in combination with --s-maxn.")
        parser.add_argument("--s-maxn", type=int, default=7100000, metavar='<int>',
                            help="Specify the number of suboptimal structures. The corresponding" +
                            " energy range is computed form the full length molecule.")

    # *barriers* arguments using argparse
    if barriers:
        parser.add_argument("--barriers", default='barriers', action='store',
                            metavar='<str>',
                            help="Specify path to your *barriers* executable")
        parser.add_argument("--b-minh", type=float, default=0.001, metavar='<flt>',
                            help="Set the minimum barrier height (i.e. barriers --minh)")
        parser.add_argument("--b-maxn", type=int, default=100, metavar='<int>',
                            help="Set the maximum number of local minima (i.e. barriers --max)")

    # *treekin* arguments using argparse
    if treekin:
        parser.add_argument("--treekin", default='treekin', action='store',
                            metavar='<str>',
                            help="Specify path to your *treekin* executable")
        # TODO: (--t0) set an appropriate default value, or None!
        parser.add_argument("--t0", type=float, default=0, metavar='<flt>',
                            help="First time point of the printed time-course")
        parser.add_argument("--ti", type=float, default=1.2, metavar='<flt>',
                            help="Output-time increment of solver (t1 * ti = t2)")
        parser.add_argument("--t8", type=float, default=0.02, metavar='<flt>',
                            help="Transcription speed in seconds per nucleotide")
        parser.add_argument("--p0", nargs='+', default=['1=1'], metavar='<int>=<flt>',
                            help="Initial population vector as a space-separated list of " +
                            "assigments \"index=occupancy\"")

    # ViennaRNA model arguments:
    if noLP:
        parser.add_argument("--noLP", action="store_true",
                            help="The ViennaRNA --noLP option")
    if temperature:
        parser.add_argument("-T", "--temperature", type=float, default=37.0,
                            metavar='<flt>', help="The temperature for ViennaRNA computations")
    if circ:
        parser.add_argument("--circ", action="store_true",
                            help="Circular RNA")

    # Other, commonly used arguments:
    if tmpdir:
        parser.add_argument("--tmpdir", default='', action='store',
                            metavar='<str>',
                            help="Specify path for storing temporary output files. " +
                            "These files will not be removed when the program terminates.")
    if name:
        parser.add_argument("--name", default='', metavar='<str>',
                            help="Name your output files, this option overwrites the fasta-header")
    if force:
        parser.add_argument("-f", "--force", action="store_true",
                            help="Force to overwrite existing BarMap files.")
    if verbose:
        parser.add_argument("-v", "--verbose", action='count', default=0,
                            help="Track process by writing verbose output during calculations.")

    # Convenient arguments for cotranscriptional folding
    if start:
        parser.add_argument("--start", type=int, default=1, metavar='<int>',
                            help="Start transcription at this nucleotide")
    if stop:
        parser.add_argument("--stop", type=int, default=None, metavar='<int>',
                            help="Stop transcription at this nucleotide")

    if k0:
        # TODO: barriers does not support this argument, it is a hack after all.
        # there is a --times argument for treekin which can be used instead.
        parser.add_argument("--k0", type=float, default=2e5, metavar='<flt>',
                            help="Arrhenius rate prefactor")

    if tX:
        parser.add_argument("--tX", type=float, default=60, metavar='<flt>',
                            help="Simulation time after transcription")

    if cutoff:
        parser.add_argument("--occupancy-cutoff", type=float, default=0.01, metavar='<flt>',
                            help="Occupancy cutoff for population transfer.")

def make_pair_table(ss, base=0, chars=['.']):
    """ Return a secondary struture in form of pair table.

    Args:
      ss (str): secondary structure in dot-bracket format
      base (int, optional): choose between a pair-table with base 0 or 1
      chars (list, optional): a list of characters to be are ignored, default:
        ['.']

    **Example:**
       base=0: ((..)). => [5,4,-1,-1,1,0,-1]
        i.e. start counting from 0, unpaired = -1
       base=1: ((..)). => [7,6,5,0,0,2,1,0]
        i.e. start counting from 1, unpaired = 0, pt[0]=len(ss)

    Returns:
      [list]: A pair-table
    """
    stack = []
    if base == 0:
        pt = [-1] * len(ss)
    elif base == 1:
        pt = [0] * (len(ss) + base)
        pt[0] = len(ss)
    else:
        raise RiboUtilsError(f"unexpected value in make_pair_table: (base = {base})")

    for i, char in enumerate(ss, base):
        if (char == '('):
            stack.append(i)
        elif (char == ')'):
            try:
                j = stack.pop()
            except IndexError as e:
                raise RiboUtilsError("Too many closing brackets in secondary structure")
            pt[i] = j
            pt[j] = i
        elif (char not in set(chars)):
            raise RiboUtilsError(f"unexpected character in sequence: '{char}'")
    if stack != []:
        raise RiboUtilsError("Too many opening brackets in secondary structure")
    return pt

def make_loop_index(ss):
    """Returns a list of loop indices to see where new base-pairs can be formed.
    """
    loop, stack = [], []
    nl, L = 0, 0
    for i, char in enumerate(ss):
        if char == '(':
            nl += 1
            L = nl
            stack.append(i)
        loop.append(L)
        if char == ')':
            _ = stack.pop()
            try:
                L = loop[stack[-1]]
            except IndexError:
                L = 0
    return loop

def parse_vienna_stdin(stdin, chars = 'ACUGN&', skip = '-'):
    """Parse name and sequence information from file with fasta format.

    Only one input-sequence is allowed at a time.

    Args:
      stdin (list): Input to parse, ususally :obj:`sys.stdin`
      chars (string, optional): Allowed characters in a sequence.

    Returns:
      str, str: name and sequence.
    """
    name = 'NoName'
    seq = ''
    for line in stdin:
        if re.match('>', line):
            if name != 'NoName':
                raise NotImplementedError(
                    'Only single-sequence fasta format supported!')
            else:
                name = line.strip().split()[0][1:]
        else:
            seq += line.strip()
    seq = seq.translate({ord(c): None for c in skip})
    m = re.search('[^' + chars + ']', seq)
    if m:
        raise RiboUtilsError("Does not look like RNA: ('{}' in '{}')".format(
            m.string[m.span()[0]], seq))
    return name, seq

def parse_ratefile(rfile, binary = False):
    """Return the content of a barriers rates-file.

    Args:
      rfile (str): Filename of a barriers rates file.
      binary (bool, optional): Set to True if the rates file is in binary format.
        Defaults to False.

    Returns:
      [[flt],[flt]]: A rate matrix.
    """
    if binary:
        with open(rfile, 'rb') as rf:
            dim, = unpack('i', rf.read(calcsize('i')))
            rm = []
            for e in range(dim):
                col = []
                for e in range(dim):
                    r, = unpack('d', rf.read(8))
                    col.append(r)
                rm.append(col)
            RM = list(map(list, list(zip(*rm))))
    else:
        RM = []
        with open(rfile) as rates:
            for line in rates:
                RM.append((list(map(float, line.strip().split()))))
    return RM

def networkx_graph(RL):
    try:
        import networkx as nx
    except ImportError as err:
        raise SystemExit(f'need networkx installed: {err}')
    pass
