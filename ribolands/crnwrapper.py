#
#  Coded by: Stefan Badelt <stef@tbi.univie.ac.at>
#  University of Vienna, Department of Theoretical Chemistry
#
#  -*- Style -*-
#  Use double quotes or '#' for comments, such that single quotes are available
#  for uncommenting large parts during testing
#
#  *) do not exceed 80 characters per line
#

from __future__ import division, print_function

from builtins import str
import subprocess as sub
from ribolands.syswraps import SubprocessError
from crnsimulator import ReactionGraph


def DiGraphSimulator(CG, fname, nlist, p0, t0, t8,
                     t_lin=300,
                     t_log=None,
                     jacobian=True,  # slower, but supposably useful.
                     # force=False, # <= not implemented
                     verb=False):
    """A wrapper function for the python module: crnsimulator

    Args:
      CG (:obj:`networkx.DiGraph`): A networkx conformation graph, where nodes
        are secondary structures and edges have weights in form of transition rates.
      fname (str): Name of the file.
      nlist (list): A list of nodes in the system.
      p0 (list): vector of initial concentrations.
      t0 (float): start of simulation.
      t8 (float): end of simulation.
      t_lin (int, optional) : evenly return output on a lin-scale from t0 to t8 (*t_lin* times)
      t_log (int, optional) : evenly return output on a log-scale from t0 to t8 (*t_log* times)
      jacobian (bool, optional): Calculate the Jacobi-Matrix for differentiation.
        Not recommended, as it slows down the computations quite a bit.
      verb (bool, optional): verbose information. Defaults to False.

    Raises:
      SubprocessError: Subprocess failed with returncode: ...

    Returns:
      [str]: The name of a treekin-like nxy file.
    """

    tfile = fname + '.tkn'
    xfile = fname + '.py'

    crn = []
    for e in CG.edges():
        if not CG.node[e[0]]['active'] or not CG.node[e[1]]['active']:
            continue
        reactant = 'id_' + str(CG.node[e[0]]['identity'])
        product = 'id_' + str(CG.node[e[1]]['identity'])
        rate = CG[e[0]][e[1]]['weight']
        crn.append([[reactant], [product], rate])

    RG = ReactionGraph(crn)
    svars = ['id_' + str(x[1]['identity']) for x in nlist]

    filename, _ = RG.write_ODE_lib(sorted_vars=svars, filename=xfile)

    # Set/Adjust Parameters
    syscall = ['python', filename]
    syscall.extend(['--nxy'])
    syscall.extend(['--t0', str(t0)])
    syscall.extend(['--t8', str(t8)])
    if t_lin:
        syscall.extend(['--t-lin', str(t_lin)])
    if t_log:
        syscall.extend(['--t-log', str(t_log)])
    syscall.extend(['--p0'])
    syscall.extend(p0)
    if verb:
        print("# {} > {}".format(' '.join(syscall), tfile))

    # Do the simulation (catch treekin errors)
    with open(tfile, 'w') as tf:
        proc = sub.Popen(syscall, stdout=tf, stderr=None)
        proc.communicate(None)
        if proc.returncode:
            call = "# {} > {}".format(' '.join(syscall), tfile)
            raise SubprocessError(proc.returncode, call)
        tf.write('# of iterations: ' + str(t_lin if t_lin else t_log) + '\n')

    return tfile
