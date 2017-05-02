
import subprocess as sub

from ribolands.syswraps import SubprocessError

def module_exists(module_name):
  try:
    __import__(module_name)
  except ImportError:
    return False
  else:
    return True

def DiGraphSimulator(CG, fname, nlist, p0, t0, t8, 
    t_lin = 300,
    t_log = None,
    jacobian=True, # slower, but supposably useful.
    #force=False, # <= not implemented
    verb = False):
  """A wrapper function for the python module: crnsimulator
  Args:
    CG <>:  A "networkx" conformation graph, where nodes are secondary structures and 
            edges are transition rates.
    fname <str>: Name of the file
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

  tfile = fname+'.tkn'
  xfile = fname+'.py'

  if module_exists('crnsimulator'):
    from crnsimulator import ReactionGraph
    from crnsimulator.odelib_template import add_integrator_args
  else :
    raise RuntimeError('Need to install Module: "crnsimulator"', 
        'Download from https://github.com/bad-ants-fleet/crnsimulator')

  crn = []
  for e in CG.edges_iter() :
    if not CG.node[e[0]]['active'] or not CG.node[e[1]]['active'] :
      continue
    reactant = 'id_' + str(CG.node[e[0]]['identity'])
    product = 'id_' + str(CG.node[e[1]]['identity'])
    rate = CG[e[0]][e[1]]['weight']
    crn.append([[reactant],[product],rate])

  RG = ReactionGraph(crn)
  svars = map(lambda x: 'id_'+str(x[1]['identity']), nlist)

  filename, _ = RG.write_ODE_lib(sorted_vars=svars, filename=xfile)


  # Set/Adjust Parameters
  syscall = ['python', filename]
  syscall.extend(['--nxy'])
  syscall.extend(['--t0', str(t0)])
  syscall.extend(['--t8', str(t8)])
  if t_lin :
    syscall.extend(['--t-lin', str(t_lin)])
  if t_log:
    syscall.extend(['--t-log', str(t_log)])
  syscall.extend(['--p0'])
  syscall.extend(p0)
  if verb :
    print "# {} > {}".format(' '.join(syscall), tfile)

  # Do the simulation (catch treekin errors)
  with open(tfile, 'w') as tf:
    proc = sub.Popen(syscall,stdout=tf,stderr=None)
    proc.communicate(None)
    if proc.returncode :
      call = "# {} > {}".format(' '.join(syscall), tfile)
      raise SubprocessError(proc.returncode, call)
    tf.write('# of iterations: ' + str(t_lin if t_lin else t_log) + '\n')

  return tfile

