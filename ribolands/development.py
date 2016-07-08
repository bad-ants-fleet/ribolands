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

#roomt = vrna_model(temp=25.0)
#frank = ribo_molecule('NoName', 'AUCG', model=roomt)
#

#success = frank.sys_suboptimals()
#bfile = frank.barriers().get_barfile()
#bfile = frank.barriers().get_ratefile()

class vrna_model(object):
  """Very much like the standard ViennaRNA model details:
  defines the energy and environmental parameters"""

  # Immutable parameters
  DEFAULT_ION_CONCENTRATIONS=[{Na:1, Cl:0, Mg:0}]
  CONST_R = constants.R / constants.calorie_th / 1000 # kcal mol^-1 K^-1

  def __init__(self, 
      params = 'turner2004.par',
      dangle = 2, # [0: none, 1: ?, 2: some, 3: mathews]
      noGU = False,
      noLP = False,
      tetra = True,
      tempC = 37.0):

    # energy model parameters
    self._params = params
    self._dangle = dangle
    self._noGU = noGU
    self._noLP = noLP
    self._tetra = tetra

    # environment and physical constants
    self._tempC = tempC # Celsius
    self._tempK = 273.15 + float(tempC) # Kelvin
    self._RT = CONST_R * self._tempK # kcal/mol

    # other
    self.theophylline = False
    self.tetracycline = False

  @property
  def params(self):
    return self._params

  @property
  def dangles(self):
    return self._dangle

# class kinetic_model(vrna_model):
#   def __init__(self, model='arrhenius'):
#     if model == 'arrhenius' :
#       self.rate = arrhenius_rate
# 
#   def arrhenius_rate(e1, e2):
#     dG = e2 - e1
#     return exp(dG/kT

class ribo_molecule(vrna_model):
  """ Define basic properties of the molecule. 
    This is the partent class
  
  """
  def __init__(self, name, sequence, circ=False, model=None, ):
    self.sequence = sequence
    self.circ = circ
    if model :
      self.model = model
    else :
      self.model = super(ribolands_molecule, self).__init__(**kwargs)

  def __add__(self, other):
    pass

class vrna_suboptimals(vrna_model):
  """ A class to compute Wuchty suboptimal structures """
  def __init__(self, 
      RNAsubopt='RNAsubopt', 
      ener=None, 
      noss = 5000, 
      sort=False
      model=None):
    self._RNAsubopt = RNAsubopt
    self._ener = ener
    self._noss = noss # number of secondary structures
    self._sort = sort
    if model is None :
      super(vrna_suboptimals, self).__init__(**kwargs)
    else :

    # Data handling
    self._dict = [{}]
    self._gzip = True
    self._file = filename

  def sys_subopt_range(self):
    # call RNAsubopt iteratively to calculate an energy range
    pass

  def sys_RNAsubopt(self):
    # use subprocess to call RNAsubopt and write to self._file
    pass

class vrna_barriers(vrna_model):
  """ A class to compute barriers energy landscapes """
  def __init__(self):
    self.minh = minh
    self.maxn = maxn
    self.rates = True
    vrna_suboptimals.__init__(self, kwargs)

    self._sfile = sfile
    self._rfile = name + '.rts'
    self._bfile = name + '.bar'
    self._pfile = name + '.ps'

  def sys_barriers(self):
    pass

class vrna_kinetics(vrna_model):
  """ A class to compute kinetics on coarse-grained landscapes """
  def __init__(self):
    self.t0 = 1e-6
    self.ti = 1.02
    self.t8 = 1e10

  def sys_treekin(self):
    pass

  def numpy_mexp(self):
    pass

  def sys_sundials(self):
    pass

  def numpy_odes(self, nxgraph, sorting=None):
    pass


""" Thinktank for a new and better interface:

  ******** Example *********

  # Standard Parameter sets:
  context = ribolands_environment(...)
  model   = ribolands_energy_model(...)

  romeo = ribolands_nucleic_acid(seq, circ=False, env=context, emodel=model)


  simulation = ril.DrTransfomer(romeo)
  simulation = ril.BarMap(romeo)
  simulation = ril.cofold_barriers(romeo, julia)

  
  context = pyrna_environment(tempC=37.0, NaCL=1.0, theo=False, ...)
  model   = pyrna_energy_model(params='turner', dangle=2, noGU=False, tetraloop=True)
  
  romeo = pyrna_nucleic_acid(seq, circ=False, environment=context, emodel=model)
  
  romeo.set_environment.tempC(37.0)
  romeo.set_environment.kT(0.61)
  romeo.set_energy_model.params('turner1999')
  
  [ss, en] = romeo.mfe(tempC=15, circ=True) # <-- Sets the variable for this run only
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



