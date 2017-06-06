#
#  Coded by: Stefan Badelt <stef@tbi.univie.ac.at>
#  University of Vienna, Department of Theoretical Chemistry
#
#  -*- Style -*- 
#  Use double quotes or '#' for comments, such that single quotes are available
#  for uncommenting large parts during testing
#
#  *) do not exceed 80 characters per line
#  *) indents: 2x whitespace, no tab characters!
#
#  -*- VIM config -*- 
#  set textwidth=80
#  set ts=2 et sw=2 sts=2
#
#  -*- Content -*-
#  *) parsers for stdin, barfiles and rate-matrix

import re
import sys
from struct import pack, unpack
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

class ProgressBar(object):
  def __init__(self, clockmax):
    if clockmax < 1:
      raise Exception("Imporper input.")

    self.clock = 0
    self.clockmax = clockmax
    print('# |Progress' + ' '*(clockmax-8) + '|')
    sys.stdout.write('# |')
    
  def inc(self):
    self.clock +=1 
    if self.clock % 10 == 0:
      sys.stdout.write(str(self.clock/10))
    elif self.clock % 5 == 0 :
      sys.stdout.write(",")
    else :
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
  if RNAsubopt :
    parser.add_argument("--RNAsubopt", default='RNAsubopt', action = 'store',
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
  if barriers :
    parser.add_argument("--barriers", default='barriers', action = 'store',
        metavar='<str>',
        help="Specify path to your *barriers* executable") 
    parser.add_argument("--b-minh", type=float, default=0.001, metavar='<flt>',
        help="Set the minimum barrier height (i.e. barriers --minh)")
    parser.add_argument("--b-maxn", type=int, default=100, metavar='<int>',
        help="Set the maximum number of local minima (i.e. barriers --max)")

  # *treekin* arguments using argparse
  if treekin :
    parser.add_argument("--treekin", default='treekin', action = 'store',
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
        help="Initial population vector as a space-separated list of "+\
            "assigments \"index=occupancy\"")

  # ViennaRNA model arguments:
  if noLP :
    parser.add_argument("--noLP", action="store_true",
        help="The ViennaRNA --noLP option")
  if temperature :
    parser.add_argument("-T","--temperature", type=float, default=37.0,
        metavar='<flt>', help="The temperature for ViennaRNA computations")
  if circ :
    parser.add_argument("--circ", action="store_true",
     help="Circular RNA")  

  # Other, commonly used arguments:
  if tmpdir :
    parser.add_argument("--tmpdir", default='', action = 'store', 
        metavar='<str>',
        help="Specify path for storing temporary output files. " + \
            "These files will not be removed when the program terminates.")
  if name :
    parser.add_argument("--name", default='', metavar='<str>',
        help="Name your output files, this option overwrites the fasta-header")
  if force :
    parser.add_argument("-f","--force", action="store_true",
        help="Force to overwrite existing BarMap files.") 
  if verbose :
    parser.add_argument("-v","--verbose", action='count', default=0,
        help="Track process by writing verbose output during calculations.") 

  # Convenient arguments for cotranscriptional folding
  if start :
    parser.add_argument("--start", type=int, default=1, metavar='<int>',
        help="Start transcription at this nucleotide")
  if stop :
    parser.add_argument("--stop", type=int, default=None, metavar='<int>',
        help="Stop transcription at this nucleotide")

  if k0 :
    # TODO: barriers does not support this argument, it is a hack after all.
    # there is a --times argument for treekin which can be used instead.
    parser.add_argument("--k0", type=float, default=2e5, metavar='<flt>',
        help="Arrhenius rate prefactor")

  if tX :
    parser.add_argument("--tX", type=float, default=60, metavar='<flt>',
        help="Simulation time after transcription")

  if cutoff :
    parser.add_argument("--occupancy-cutoff", type=float, default=0.01, metavar='<flt>',
        help="Occupancy cutoff for population transfer.")


def make_pair_table(ss, base=0, chars=['.']):
  """Return a secondary struture in form of pair table.

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
  stack=[];

  if base is 0:
    pt=[-1] * len(ss);
  elif base == 1:
    pt = [0] * (len(ss) + base);
    pt[0] = len(ss);
  else :
    raise ValueError("unexpected value in make_pair_table: \
        (base = " + str(base) + ")")

  for i, char in enumerate(ss, base):
    if (char == '('):
      stack.append(i);
    elif (char == ')'):
      try :
        j=stack.pop();
      except IndexError, e :
        raise RuntimeError("Too many closing brackets in secondary structure")
      pt[i]=j
      pt[j]=i
    elif (char not in set(chars)):
      raise ValueError("unexpected character in sequence: '" + char + "'")

  if stack != [] :
    raise RuntimeError("Too many opening brackets in secondary structure")
  return pt

def parse_vienna_stdin(stdin, chars=['A','C','U','G','&']):
  """Parse name and sequence information from file with fasta format.
  
  Only one Input-Sequence is allowed at a time. 

  Args:
    stdin (list): Input to parse, ususally :obj:`sys.stdin`
    chars (list, optional): Allowed characters in a sequence.

  Returns:
    [(str, str)]: A tuple containing name and sequence.
  """
  name = 'NoName'
  seq  = ''
  for line in stdin:
    if re.match('>', line):
      if name != 'NoName' :
        raise NotImplementedError('Only single-sequence fasta format supported!')
      else : 
        name = line.strip().split()[0][1:]
    else:
      seq += line.strip()

  m = re.search('[^'+''.join(chars)+']', seq) 
  if m :
    raise ValueError("Does not look like RNA: ('{}' in '{}')".format(m.string[m.span()[0]], seq))
  return (name, seq)

def parse_barfile(bfile, seq=''):
  """Return the content of a barriers output-file. 

  Args:
    bfile (str): Filename of a barriers output file.
    seq (str, optional): Raises an error if the sequence is provided and
      different than the one in the barfile.
  
  Raises:
    ValueError: Wrong sequence.

  Returns:
    [list of lists]: The content of the barfile.
  """
  # NOTE: This function would make more sense if it returns an object that
  # stores barriers info in a consistent way.
  output = []
  with open(bfile) as bar :
    for e, line in enumerate(bar) :
      if e == 0 : 
        if seq and seq != line.strip() :
          raise ValueError('Wrong sequence ' + seq + ' vs. ' + line)
      else :
        output.append(line.strip().split())
        #[idx, lmin, en, father, bar] = line.strip().split()[0:5]
        #output.append([idx, lmin, en, father, bar])
  return output

def parse_ratefile(rfile, binary=False):
  """Return the content of a barriers rates-file. 

  Args:
    rfile (str): Filename of a barriers rates file.
    binary (bool, optional): Set to True if the rates file is in binary format.
      Defaults to False.
  
  Returns:
    [[flt],[flt]]: A rate matrix.
  """
  if binary :
    with open(rfile) as rf:
      dim, = unpack('i', rf.read(4))
      rm = []
      for e in range(dim):
        col = []
        for e in range(dim):
          r, = unpack('d', rf.read(8))
          col.append(r)
        rm.append(col)
      RM = map(list, zip(*rm))
  else :
    RM = []
    with open(rfile) as rates :
      for line in rates :
        RM.append((map(float, line.strip().split())))

  return RM

def plot_nxy(name, tfile, 
    title='',
    plim=1e-2,
    lines=[],
    ylim=(None,None),
    xlim=(None,None), 
    lilog=None):
  """ Plot a list of trajectories.

  Args:
    name (str): Name of the pdf file.
    tfile (str): Filename of a treekin `nxy` output file.
    title (str, optional): Name of the title for the plot.
    plim (float, optional): Minimal occupancy to plot a trajectory. Defaults to 0.01
    lines ([int,..], optional): Selected list of lines to be plotted.
    ylim ((float,float), optional): matplotlib ylim.
    xlim ((float,float), optional): matplotlib xlim.
    lilog (float, optional): Set a divider of the x-axis into a linear and
      logarithmic part. Requires xlim and ylim to be specified.

  Returns:
    [str]: Name of the output file.
  """
  lines=set(lines)
  title = title if title else name

  nxy=[]
  with open(tfile) as tkn :
    for line in tkn :
      if re.match('#', line): 
        continue
      nxy.append(map(float, line.strip().split()))

  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  if None not in ylim : ax.set_ylim(ylim)
  if None not in xlim : ax.set_xlim(xlim) 

  if lilog :
    if None in xlim : 
      raise ValueError('Specify xlim')
    if None in ylim : 
      raise ValueError('Specify ylim')

    # Make the second part of the plot logarithmic
    ax.set_xscale('linear')
    ax.set_xlim((xlim[0], lilog))
    divider = make_axes_locatable(ax)

    axLog = divider.append_axes("right", size=2, pad=0.0, sharey=ax)
    axLog.set_xscale('log')
    axLog.set_xlim((lilog+0.00001, xlim[1]))
    axLog.set_ylim(ylim)
  else :
    ax.set_xscale('log')

  time = []
  for e, traject in enumerate(zip(*nxy)):
    if e==0: 
      time = traject
      continue
    if lines and e not in lines : continue
    if plim and max(traject) < plim : continue
    if ylim and max(traject) > ylim : continue
    p, = ax.plot(time, traject, '-')
    if lilog: p, = axLog.plot(time, traject, '-')
    p.set_label("lmin {:d}".format(e))
    #p.set_color(axcolors[finalmin])

  fig.set_size_inches(7,3)
  fig.text(0.5,0.95, title, ha='center', va='center')
  plt.xlabel('time [seconds]', fontsize=11)
  plt.ylabel('occupancy [mol/l]', fontsize=11)

  # """ Add ticks for 1 minute, 1 hour, 1 day, 1 year """
  # plt.axvline(x=60, linewidth=1, color='black', linestyle='--')
  # plt.axvline(x=3600, linewidth=1, color='black', linestyle='--')
  # plt.axvline(x=86400, linewidth=1, color='black', linestyle='--')
  # plt.axvline(x=31536000, linewidth=1, color='black', linestyle='--')
  plt.legend()

  plt.savefig(name, bbox_inches='tight')
  return name

