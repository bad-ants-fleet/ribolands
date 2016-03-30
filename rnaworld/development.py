#!/usr/bin/env python

import re
import sys
import os
import argparse
import collections as c
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# call: 
#   odes = ode_system(CG, [options, sorting, ...]) 
#   jaco = jacobian(odes)
#   solver = sys_sundials(ode, jac=jaco, kwargs)
#   python_solver(ode, jac=jaco, kwargs, plot=True)


def ode_solver(GraphX, 
    variables=None, var_vert=None, vert_var=None, verb=True):

  return 


""" Thinktank for a new and better interface:

  ******** Example *********
  
  context = pyrna_environment(tempC=37.0, NaCL=1.0, theo=False, ...)
  model   = pyrna_energy_model(params='turner', dangle=2, noGU=False, tetraloop=True)
  
  romeo = pyrna_nucleic_acid(seq, circ=False, environment=context, emodel=model)
  
  romeo.set_environment.tempC(37.0)
  romeo.set_environment.kT(0.61)
  romeo.set_energy_model.params('turner1999')
  
  [ss, en] = romeo.mfe(tempC=15, circ=True) # Sets the variable for this run only
  # alternative: pyrna_mfe(romeo, tempC=15, circ=True)
  efe = romeo.efe
  pfc = romeo.pfc
  ffile = romeo.sys_fold
  
  romeo.energy_model.noLP(True)
  romeo.energy_model.k0(2e5)
  sfile = romeo.sys_suboptimals(fname=romeo, ener=None, maxn=5000, maxe=30.00)
  [bfile, rfile, ... ] = romeo.sys_barriers(minh=1, maxn=50, rates=True, k0=1)
  tfile = romeo.sys_treekin(t0=1e-6, ti=1.02, t8=1e6)
  
  what about two sequences?
  context = romeo.get_environment()
  model = romeo.get_energy_model()
  julia = pyrna_nucleic_acid(seq, context=context, model=model)
  
  
  couple = romeo + julia # only works if context is the same (?)
  
  cut_point = couple.get_cut_points
  couple.mfe(cut_point=-1)
  copule.sys_suboptimals
  
  onemore += felix # will be prohibited, unless using NUPACK

  ******** Structure *********

  Initialize parent objects:
    * vrna_environment # temp, circ
    * energy_model # 
    * model_options(subopt, barriers, treekin)
  
  Initialize RNA object and inherit from parents, if they are set:
    * pyrna_nucleic_acid(seq, env=env, emod=emod)
  
  Call functions with obligatory one-time parameters

"""

""" under construction """
class pyrna_context:
  """ Set model details to standard vrna details """
  temperature = 37.0
  theophylline = None
  tetracycline = None
  ions = None

# inhertiance: children have access to parent methods
class pyrna_molecule: 
  """ Define basic properties of the molecule. 
    This is the partent class
  
  """
  def __init__(self, sequence):
    self.sequence = sequence
    circ = False
    print "Creating an RNA object"

  def energy_model(self):
    parameter_file = 'turner1999'
    dangles = 2 # [0: none, 1: ?, 2: some, 3: mathews]
    GUclosing = True
    tetraloopen = True


class pyrna_options:
  def __init__(self):
    pass

  def general(self):
    noLP = False

  def suboptimals(self):
    subopt_energy = None
    maxspt_energy = 30.0
    subopt_structs= 5000

  def barriers(self):
    rates=True
    bsize=False
    saddle=False
    minh = 0.01
    maxn = 100

  def kinetics(self):
    t0 = 1e-6
    ti = 1.02
    t8 = 1e10


""" Under Testing """

class cd:
  """ Context manager for changing the current working directory 
    http://stackoverflow.com/questions/431684/how-do-i-cd-in-python 
  """
  def __init__(self, newPath):
    self.newPath = os.path.expanduser(newPath)

  def __enter__(self):
    self.savedPath = os.getcwd()
    os.chdir(self.newPath)

  def __exit__(self, etype, value, traceback):
    os.chdir(self.savedPath)

'''
with cd("~/Library"):
  s.call("ls")
'''

def model_details(force=False, tmp=27, **kwars):
  def wrapper(name):
    return "bla-{0}".format(func(name))
  return wrapper

""" utilities """

def barriersCG(mfile, efile, verb=False):
  """ Takes a sequence and a file with structures, populations, ...
  (mfile). Computes a barrier tree for the given sequence and returns 
  a dictionary with the correctly mapped populations:
  (requires barriers-v1.6 or later)

  mfile in form:
  structure energy population (+ lost columns)
  ...
  
  :return: Dictionary[lmin] = population
  """
  
  BarCG = c.defaultdict(float)
  with open(mfile, 'r') as m, open(efile, 'r') as mapstruc:
    for line in mapstruc:
      if re.match('[\.\(\)]+', line.strip()):
        gstr, sptidx, energy, fmin, fminT, gmin, gminT = line.strip().split()
        [old, pop, en] = m.readline().strip().split()[0:3]
        if verb : print [gminT, gstr, old, en, pop]
        BarCG[int(gminT)] += float(pop)
      elif re.match('not in hash', line.strip()):
        [old, pop, en] = m.readline().strip().split()[0:3]
        print >> sys.stderr, old, en, "structure not in hash"
        sys.exit('over and out')
      elif re.match('not yet assigned', line.strip()):
        [old, pop, en] = m.readline().strip().split()[0:3]
        print >> sys.stderr, old, en, "structure not yet assigned"
        sys.exit('over and out')
  return BarCG



