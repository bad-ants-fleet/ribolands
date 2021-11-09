#
# ribolands.crnwrapper
#
# TODO: 
#   - needs testing.
#   - add force option?
#   - make a more general interface that works with all types landscape objects?
#

import logging
rlog = logging.getLogger(__name__)

import subprocess as sub
from ribolands.syswraps import SubprocessError
from crnsimulator import ReactionGraph

def DiGraphSimulator(CG, oname, nlist, p0, t0, t8,
                     t_lin = 300,
                     t_log = None):
    """ A wrapper function for the python module: crnsimulator

    Args:
      CG (:obj:`networkx.DiGraph`): A networkx conformation graph, where nodes
        are secondary structures and edges have weights in form of transition rates.
      oname (str): Name the output file.
      nlist (list): A list of nodes in the system.
      p0 (list): vector of initial concentrations.
      t0 (float): start of simulation.
      t8 (float): end of simulation.
      t_lin (int, optional) : evenly return output on a lin-scale from t0 to t8 (*t_lin* times)
      t_log (int, optional) : evenly return output on a log-scale from t0 to t8 (*t_log* times)

    Raises:
      SubprocessError: Subprocess failed with returncode: ...

    Returns:
      [str]: The name of a treekin-like nxy file.
    """

    tfile = oname + '.tkn'
    xfile = oname + '.py'

    crn = []
    for e in CG.edges():
        if not CG.nodes[e[0]]['active'] or not CG.nodes[e[1]]['active']:
            continue
        reactant = 'id_' + str(CG.nodes[e[0]]['identity'])
        product = 'id_' + str(CG.nodes[e[1]]['identity'])
        rate = CG[e[0]][e[1]]['weight']
        crn.append([[reactant], [product], rate])

    RG = ReactionGraph(crn)
    svars = ['id_' + str(CG.nodes[x]['identity']) for x in nlist]

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
    rlog.debug(f"{' '.join(syscall)} > {tfile}")

    # Do the simulation (catch treekin errors)
    with open(tfile, 'w') as tf:
        proc = sub.Popen(syscall, stdout=tf, stderr=None)
        proc.communicate(None)
        if proc.returncode:
            call = "# {} > {}".format(' '.join(syscall), tfile)
            raise SubprocessError(proc.returncode, call)
        tf.write('# of iterations: ' + str(t_lin if t_lin else t_log) + '\n')

    return tfile
