
import imp
import sympy
import numpy as np
import networkx as nx # Not needed here, but it is the Graph format used
from scipy.integrate import odeint

from crnsimulator import writeODElib
from crnsimulator.reactiongraph import DiGraph_to_ODE

def DiGraphSimulator(CG, _fname, _tmpdir, _odename, nlist, p0, t0, t8, 
    t_lin = 300,
    t_log = None,
    jacobian=True, # slower, but supposably useful.
    force=False, 
    verb = False):
  """A wrapper function for the python module: crnsimulator
  Args:
    CG <>:  A "networkx" conformation graph, where nodes are secondary structures and 
            edges are transition rates.
    _fname <str>: Name of the file
    _tmpdir <str>: Directory to store temporary data files.
    _odename <str>: Name of the ODE system
    nlist <lol>: A list of nodes in the system (redundant?)
    p0 <list>: vector of initial concentrations
    t0 <float>: start of simulation
    t8 <float>: end of simulation
    t_lin <int> : evenly return output on a lin-scale from t0 to t8 (*t_lin* times)
    t_log <int> : evenly return output on a log-scale from t0 to t8 (*t_log* times)
    jacobian <bool>: Calculate the Jacobi-Matrix for differentiation. Not recommended, 
                      as it slows down the computations quite a bit.
    verb <bool>: print a commandline call to run the simulation (debugging)
    
  """

  crn, ode, rdict = DiGraph_to_ODE(CG, rate_dict=True)

  M = []
  V = []
  for (ss, attr) in nlist:
    ssid = 'id_' + str(attr['identity'])
    sfunc = sympy.sympify(' + '.join(
      ['*'.join(map(str,xp)) for xp in ode[ssid]]))
    ode[ssid] = sfunc
    M.append(sfunc)
    V.append(ssid)

  M = sympy.Matrix(M)

  if jacobian:
    # NOTE: the sympy version breaks regularly:
    # J = M.jacobian(V) 
    # This is why it's done per pedes:
    J = []
    for f in M :
      for x in V:
        J.append(f.diff(x))
    J = sympy.Matrix(J)
  else :
    J = None

  ### => SOLVER
  odefile =  writeODElib(V, M, odename=_odename, path=_tmpdir, jacobian=J, rdict=rdict)
  if verb :
    print '# Wrote ODE system:', odefile

  _temp = imp.load_source(_odename, odefile)
  odesystem = getattr(_temp, _odename)
  if J :
    print 'using jacobian!'
    jacobian = getattr(_temp, 'jacobian')
  else :
    jacobian = None

  if t_lin:
    time = np.linspace(t0, t8, t_lin)
  elif t_log:
    time = np.logspace(t0, t8, t_log)
  else :
    raise ValueError('Need to set either t_lin or t_log!')

  myp0 = [0] * len(nlist)
  for tup in p0 :
    p, f = tup.split('=')
    myp0[int(p)-1]=float(f)

  # Set/Adjust Parameters
  tfile = _fname+'.tkn'
  rates = None # use the rates from odefile
  try :
    if verb :
      syscall = ['python', odefile]
      syscall.extend(['--nxy'])
      syscall.extend(['--noplot'])
      syscall.extend(['--t0', str(t0)])
      syscall.extend(['--t8', str(t8)])
      syscall.extend(['--t-lin', str(t_lin if t_lin else t_log)])
      syscall.extend(['--p0'])
      syscall.extend(p0)
      print "# {} > {}".format(' '.join(syscall), tfile)

    if jacobian :
      ny = odeint(odesystem, myp0, time, (rates, ), Dfun=jacobian).T
    else :
      ny = odeint(odesystem, myp0, time, (rates, )).T
  except ValueError, e:
    print odeint(odesystem, myp0, time, (rates, ), full_output=True)
    raise SystemExit('Crnwrapper failed: please submit a bug-report!')

  with open(tfile, 'w') as tf :
    for i in zip(time, *ny):
      tf.write(' '.join(map("{:.9e}".format, i))+'\n')
    tf.write('# of iterations: ' + str(t_lin if t_lin else t_log) + '\n')

  return tfile

