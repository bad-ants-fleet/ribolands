# -*- coding: utf-8 -*-

__version__ = "0.2.0"

_MIN_TREEKIN_VERSION = "0.4.1"
_MIN_BARRIERS_VERSION = "1.6.0"
_MIN_VIENNARNA_VERSION = "2.2.0"

from utils import argparse_add_arguments
from utils import parse_vienna_stdin
from utils import parse_barfile
from utils import parse_ratefile
from utils import make_pair_table
from utils import plot_simulation
from utils import ProgressBar

from syswraps import sys_subopt_range
from syswraps import sys_suboptimals 
from syswraps import sys_barriers
from syswraps import sys_treekin

from crnwrapper import DiGraphSimulator

#'''
#def getversion():
#  v_utils = utils.__version__
#  print v_utils
#'''

