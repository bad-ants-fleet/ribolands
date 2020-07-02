#
# ribolands.ribolands
# 
# Home of the RiboLandscape object.
#

import logging
rlog = logging.getLogger(__name__)

import RNA 
import math
import networkx as nx
import matplotlib.pyplot as plt
from struct import pack
from itertools import product, combinations, permutations
from crnsimulator import ReactionGraph
from crnsimulator.crn_parser import parse_crn_string

from ribolands.utils import natural_sort
from ribolands.pathfinder import apply_bp_change, get_fpath_flooding_cache

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Custom error definitions                                                     #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
class RiboLandscapeError(Exception):
    pass

class RiboLandscape(nx.DiGraph):
    """ Implemented for unimolecular reactions.

    This is the backbone of any ribolands landscape object.

    Suggested attributes for nodes: energy, structure, identity, active
    Suggested attributes for edges: weight, saddleE.
    """

    def __init__(self, sequence, vrna_md = None):
        super(RiboLandscape, self).__init__()
        self.sequence = sequence
        if vrna_md:
            self.md = vrna_md
            self.fc = RNA.fold_compound(sequence, vrna_md)
        else:
            self.md = None
            self.fc = None

        self.nodeID = 0  # A counter for autogenerated node IDs

        # Store coarse graining results.
        self.minh = None # kcal/mol; the minimal height of an energy barrier.
        self.k0 = 1
        self.lminreps = None
        self.hiddennodes = None

    def addnode(self, key, 
            structure = None, 
            identity = None, 
            energy = None, 
            occupancy = 0, 
            active = None,
            **kwargs):
        # Wrapper to add a secondary structure to the graph.
        if structure is not None and energy is None:
            energy = round(self.fc.eval_structure(structure), 2)
        if identity is None:
            identity = self.nodeID
        assert isinstance(energy, float) if energy else True
        self.add_node(key,
                identity = identity,
                structure = structure, 
                energy = energy,
                occupancy = occupancy, 
                active = active,
                **kwargs)
        self.nodeID += 1
   
    def addedges(self, n1, n2, fwb, rvb, sE):
        assert self.has_node(n1)
        assert self.has_node(n2)
        self.add_edge(n1, n2, weight = self.k0 * math.e**(-fwb/self.RT), saddleE = sE)
        self.add_edge(n2, n1, weight = self.k0 * math.e**(-rvb/self.RT), saddleE = sE)

    def edgeupdate(self, n1, n2, fwb, rvb, sE):
        self[n1][n2]['weight'] = self.k0 * math.e**(-fwb/self.RT)
        self[n2][n1]['weight'] = self.k0 * math.e**(-rvb/self.RT)
        self[n1][n2]['saddleE'] = sE
        self[n2][n1]['saddleE'] = sE

    @property
    def RT(self):
        RT = 0.61632077549999997
        if self.md.temperature != 37.0:
            kelvin = 273.15 + self.md.temperature
            RT = (RT / 310.15) * kelvin
        return RT

    @property
    def active_subgraph(self):
        """ Return a new object with only active nodes (and edges) """
        other = self.__class__(self.sequence)
        other.md = self.md
        other.fc = self.fc
        other.minh = self.minh

        for n, data in self.nodes.data():
            if data['active'] is None:
                raise RiboLandscapeError("Did you forget to call coarse-graining",
                        "before extracting the active subgraph?")
            if data['active']:
                other.add_node(n, **data)
        for n1, n2, data in self.edges.data():
            if self.nodes[n1]['active'] and self.nodes[n2]['active']:
                other.add_edge(n1, n2, **data)
        return other
    
    @property
    def active_nodes(self):
        return [n for n in self.nodes if self.nodes[n]['active'] is True]

    @property
    def inactive_nodes(self):
        return [n for n in self.nodes if self.nodes[n]['active'] is False]

    @property
    def new_nodes(self):
        return [n for n in self.nodes if self.nodes[n]['active'] is None]

    def sorted_nodes(self, attribute = 'energy', rev = False, nodes = None):
        """ Provide active nodes or new nodes, etc. if needed. """
        if nodes is None:
            nodes = self.nodes
        return sorted(nodes, key = lambda x: self.nodes[x][attribute], reverse = rev)

    def get_saddle(self, s1, s2):
        """ Returns the saddle energy of a transition edge. """
        return self[s1][s2]['saddleE'] if self.has_edge(s1, s2) else None

    def get_barrier(self, s1, s2):
        """ Returns the barrier energy of a transition edge. """
        if self.has_edge(s1, s2):
            return self[s1][s2]['saddleE'] - self.nodes[s1]['energy']
        else:
            return None

    def get_rate(self, s1, s2):
        """ Returns the direct transition rate of two secondary structures. """
        return self[s1][s2]['weight'] if self.has_edge(s1, s2) else 0

    def get_simulation_files_tkn(self, basename, snodes = None):
        """ Print a rate matrix and the initial occupancy vector.

        This function prints files and parameters to simulate dynamics using the
        commandline tool treekin. A *.bar file contains a sorted list of present
        structures, their energy and their neighborhood and the corresponding
        energy barriers. A *.rts or *.rts.bin file contains the matrix of
        transition rates either in text or binary format. Additionaly, it returns
        a vector "p0", which contains the present occupancy of structures. The
        order or elements in p0 contains

        Note:
          A *.bar file contains the energy barriers to all other nodes in the
          graph. So it is not the same as a "classic" barfile produce by
          barriers.

        Args:
            basename (str): Basename of output files.

        Returns:
            [str, str, str]: Binary rates file, Text rates file, barriers-like output.
        """
        seq = self.sequence
        if snodes is None:
            snodes = self.sorted_nodes()
        num = len(snodes) + 1

        bofile = basename + '_lands.bar'
        brfile = basename + '_rates.txt'
        bbfile = basename + '_rates.bin'
        p0 = []

        with open(bofile, 'w') as bar, open(brfile, 'w') as rts, open(bbfile, 'wb') as brts:
            bar.write("  ID {}  Energy  {}\n".format(seq, 
                ' '.join(map("{:7d}".format, range(1, num)))))
            brts.write(pack("i", len(snodes)))
            for ni, node in enumerate(snodes, 1):
                ns = self.nodes[node]['structure']
                ne = self.nodes[node]['energy']
                no = self.nodes[node]['occupancy']
                
                # Calculate barrier heights to all other basins.
                barstr = ''
                for other in snodes:
                    os = self.nodes[other]['structure']
                    oe = self.nodes[other]['energy']
                    sE = self.get_saddle(node, other)
                    if sE is not None:
                        barstr += ' {:7.2f}'.format(sE - ne)
                    else:
                        barstr += ' {:7.2f}'.format(float('nan'))

                # Print structures and neighbors to bfile:
                bar.write("{:4d} {} {:7.2f} {}\n".format(ni, ns, ne, barstr))

                # Add ni occupancy to p0
                if no > 0:
                    p0.append("{}={}".format(ni, no))

                # Print rate matrix to rfile and brfile
                trates = []
                rates = []
                for other in snodes:
                    if self.has_edge(node, other):
                        rates.append(self[node][other]['weight'])
                    else:
                        rates.append(0)

                    if self.has_edge(other, node):
                        trates.append(self[other][node]['weight'])
                    else:
                        trates.append(0)
                line = " ".join(map("{:10.4g}".format, rates))
                rts.write("{}\n".format(line))
                for r in trates:
                    brts.write(pack("d", r))
        return bbfile, brfile, bofile, p0

    def to_crn(self, filename = None, prefix = 'ID_'):

        def crn_gen():
            for (x,y) in self.edges:
                if isinstance(self.nodes[x]['identity'], int):
                    reactant = '{:s}{:d}'.format(prefix, self.nodes[x]['identity'])
                else:
                    reactant = self.nodes[x]['identity']

                if isinstance(self.nodes[y]['identity'], int):
                    product = '{:s}{:d}'.format(prefix, self.nodes[y]['identity'])
                else:
                    product = self.nodes[y]['identity']

                if self[x][y]['weight'] == 0:
                    continue
                yield "{:s} -> {:s} [k = {:g}]\n".format(reactant, product, self[x][y]['weight'])

        if filename:
            with open(filename, 'w') as crn:
                for x in crn_gen():
                    crn.write(x)
            return
        else:
            return "".join(crn_gen())

    def to_crnsimulator(self, filename, sorted_vars = None):
        """
        """
        if filename[-3:] != '.py':
            filename += '.py' 

        def irrev(rev):
            # Split CRN into irreversible reactions
            new = []
            for [r, p, k] in rev:
                if None in k:
                    print('# Set missing rates to 1.')
                    k[:] = [x if x is not None else 1 for x in k]

                if len(k) == 2:
                    new.append([r, p, k[0]])
                    new.append([p, r, k[1]])
                else:
                    new.append([r, p, k[0]])
            return new

        crnstring = self.to_crn()
        crn, species = parse_crn_string(crnstring)
        crn = irrev(crn)
        RG = ReactionGraph(crn)
        V = sorted_vars if sorted_vars else natural_sort(species)

        # ********************* #
        # PRINT ODE TO TEMPLATE #
        # ..................... #
        filename, _ = RG.write_ODE_lib(sorted_vars = V, concvect = None,
                                             jacobian = False,
                                             filename = filename)
        rlog.info(f'# Wrote ODE system: {filename}')
        return filename

    def simulate_crn(self, filename = None, sorted_vars = None):
        # Take odeint parameters
        # Take plotting parameters
        # call to_crnsimulator and do the I/O handling
        raise NotImplementedError

    def plot_to(self, oname, label = 'identity', nodes = None):
        plt.subplot(111)
        if nodes is None:
            nodes = self.nodes
        labs = {a : b[label] for a, b in self.nodes(data = True)} if label else None
        nx.draw_circular(self, with_labels = True, labels = labs, nodelist = nodes)
        plt.savefig(oname)
        plt.close()

    def coarse_grain(self, minh = None):
        """Take all (in)active nodes and coarse grain them into in/active nodes.

        Note: During that process, inactive nodes can become active and vice
        versa. Also, the internal dictionaries lminreps and hiddennodes are
        overwritten. If you want to ignore inactive nodes in the
        coarse-graining process, then extract the active graph first.
        """
        if minh is None: minh = self.minh

        self.lminreps = dict()    # what we care about ... where did a node end up?
        self.hiddennodes = dict() # bookkeeping, what nodes are part of this lmin?
        for node in self.sorted_nodes(rev = True):
            reps = self.merge_to_representatives(node, minh)
            assert self.nodes[node]['active'] is not None
            if len(reps):
                self.lminreps[node] = reps # the node was merged to reps
                for rep in reps:      # representatives are lmins (for now)!
                    if rep in self.hiddennodes:
                        self.hiddennodes[rep].add(node)
                    else:
                        self.hiddennodes[rep] = set([node])
                if node in self.hiddennodes: # if a merged node was considered an lmin earlier ...
                    basin = self.hiddennodes[node]
                    for hn in basin:
                        assert hn in self.lminreps
                        self.lminreps[hn] |= reps 
                        self.lminreps[hn] -= set([node])
                        for rep in reps:
                            self.hiddennodes[rep].add(hn)
                    del self.hiddennodes[node]
            #else: # the node was merged to itself
            #    self.lminreps[node] = set([node]) 
            #    if node not in self.hiddennodes:
            #        self.hiddennodes[node] = set([node]) # the node was merged to itself
            #    else:
            #        self.hiddennodes[node].add(node) # the node was merged to itself
        return self.lminreps, self.hiddennodes

    def merge_to_representatives(self, node, minh = None, mode = 'flooding'):
        """Identify one or more local minimum representatives for a node and
        connect all node-neighbors to those representatives.

        NOTE: mode = 'absolute' is for testing purposes only.

        NOTE: If the node finds a lmin representative with equal energy, then that
        neighbor will become active and node become inactive. This seems necessary
        because we do not want to restrict the neighbors to active conformations, 
        since that might find the wrong local minima representatives.
        """
        if minh is None:
            minh = self.minh
        assert minh is not None
        en = self.nodes[node]['energy']

        # Sort neighbors from lowest to highest energy
        # NOTE: A neighbor can be an inactive node that becomes active in the process..
        # This ensures that we always find the lowest-energy representative.
        nbrs = sorted(self.successors(node), 
                key = lambda x: (self.nodes[x]['energy'], x), reverse = False)

        # starting maximum barrier is just a tick lower than minh
        (reps, minsE) = (set(), en + minh - 0.01) 
        for nbr in nbrs:
            if self.nodes[nbr]['energy'] >= en + 0.0001:
                break
            sE = self.get_saddle(node, nbr)
            if mode == 'flooding':
                if sE + 0.01 < minsE : 
                    # a neighbor with truly lower saddle-energy 
                    (reps, minsE) = (set([nbr]), sE)
                elif abs(sE - minsE) < 0.0001 : 
                    # a neighbor with equally best saddle-energy 
                    reps.add(nbr)
            elif mode == 'absolute':
                if sE <= minsE : # a valid neighbor 
                    reps.add(nbr)

        # Connect the representatives with all other neighbors
        if len(reps):
            for lm in reps:
                self.nodes[lm]['active'] = True
                self.nodes[lm]['last_seen'] = 0
                for nbr in nbrs: 
                    if nbr == lm: 
                        continue
                    sE = max(self.get_saddle(nbr, node), self.get_saddle(node, lm)) 
                    if self.has_edge(lm, nbr) and self.get_saddle(lm, nbr) < sE:
                        continue
                    fw = sE - self.nodes[lm]['energy']
                    rv = sE - self.nodes[nbr]['energy']
                    self.addedges(lm, nbr, fw, rv, sE)
            self.nodes[node]['active'] = False
        else:
            self.nodes[node]['active'] = True
            self.nodes[node]['last_seen'] = 0
        return reps
 
class PrimePathLandscape(RiboLandscape):
    """
    An extension of the RiboLandscape object to handle prime path decomposition.
    """

    def __init__(self, *kargs, **kwargs):
        super(PrimePathLandscape, self).__init__(*kargs, **kwargs)

    def connect_nodes_n2(self, nodes = None, **kwargs):
        return self.connect_nodes_nm(nodesA = nodes, nodesB = nodes, **kwargs)

    def connect_to_active(self, nodes = None, **kwargs):
        active = self.active_nodes
        return self.connect_nodes_nm(nodesA = active, nodesB = nodes, **kwargs)

    def connect_nodes_nm(self, nodesA = None, nodesB = None, direct = False):
        """ Connect all by all, no sorting, no stopping condition.

        Args:
            nodesA (list, optional): specify a subset of nodes (e.g. new_nodes, 
                active_nodes, ...). Defaults to None: all nodes.
            nodesB (list, optional): specify a subset of nodes (e.g. new_nodes, 
                active_nodes, ...). Defaults to None: all nodes.
            direct (bool, optional): Only return paths composed of a single prime path step.

        Returns:
            list: A list of new nodes.
        """
        if nodesA is None:
            nodesA = self.nodes
        if nodesB is None:
            nodesB = self.nodes

        newnodes = set()
        for ss1, ss2 in product(nodesA, nodesB):
            if ss1 == ss2: continue
            primepath = self.get_prime_path_minima(ss1, ss2, single_step = direct)
            if primepath is None:
                assert direct is True
                continue
            for start, stop, fwb, fwg, rvb in primepath:
                assert start in self.nodes
                if stop not in self.nodes:
                    myen = round(self.nodes[start]['energy'] + fwg, 2)
                    self.addnode(stop, structure = stop, energy = myen)
                    newnodes.add(stop)
        
                sE = round(self.nodes[start]['energy'] + fwb, 2)
                assert sE == round(self.nodes[stop]['energy'] + rvb, 2)
        
                if self.has_edge(start, stop):
                    if sE < self.get_saddle(start, stop):
                        self.edgeupdate(start, stop, fwb, rvb, sE)
                else:
                    self.addedges(start, stop, fwb, rvb, sE)
        return newnodes

    def get_prime_path_minima(self, ss1, ss2, single_step = False):
        """ Use the prime-path model to find all macrostate-paths between two nodes.

        Args:
            ss1 (string): Secondary structure 1
            ss2 (string): Secondary structure 2
            single_step (bool, optional): Return the prime path only if it is just 
                a single step. Defaults to False.
        """
        assert self.minh is not None
        seq = self.sequence

        valid = [True]
        primepath = []
        def get_ppms(seq, ss1, ss2, startseq, start, stop):
            # get prime path minima
            rlog.debug(f'My input: \n{seq}\n{ss1}\n{ss2}')
            lmp = get_fpath_flooding_cache(seq, ss1, ss2, self.md, self.minh)
            rlog.debug(f'My lminpath: {lmp}')
            if lmp[0] is None:
                if isinstance(lmp[1], set):
                    backup = start
                    for x in permutations(lmp[1]):
                        start = backup
                        for (lsq, pm1, pm2) in x: 
                            get_ppms(lsq, pm1, pm2, startseq, start, stop)
                            start = apply_bp_change(startseq, start, stop, lsq, pm1, pm2)
                elif single_step: 
                    valid[0] = False
                else:
                    for (pm1, pm2) in lmp[1]: 
                        stop = apply_bp_change(startseq, start, None, seq, pm1, pm2)
                        get_ppms(seq, pm1, pm2, startseq, start, stop)
                        start = stop
            else:
                stop = apply_bp_change(startseq, start, stop, seq, ss1, ss2)
                revp = get_fpath_flooding_cache(seq, ss2, ss1, None, None)
                primepath.append([start, stop, lmp[0], lmp[1], revp[0]])
        get_ppms(seq, ss1, ss2, seq, ss1, ss2)
        return primepath if valid[0] else None

    def minimal_prime_path_graph(self, nodes = None):
        """ Uses lminreps.

        Hmmmm... so what is the difference of just taking out the interesting
        lmins and do the add_prime_path procedure from scratch and then do the
        coarse graining? Is it possible that we get better connections if we
        have a more complete graph to begin with?

        - I think it is the same effect, unless one of the additional
          conformations pulled a path minimum into a better local minimum.
          However, typically one does not know which nodes for the minimal
          prime path graph to keep until after a simulation. So you start with
          a larger landscape, then reduce it to a minimum before expansion ...
        """
        assert len(nodes) >= 2
        lminreps = dict() if self.lminreps is None else self.lminreps

        keep = set() # All nodes we are going to keep.
        for x, y in combinations(nodes, 2):
            for [start, stop, fwb, fwg, rvb] in self.get_prime_path_minima(x, y):
                basin1 = set([start]) if start not in lminreps else lminreps[start]
                basin2 = set([stop]) if stop not in lminreps else lminreps[stop]
                # NOTE: pruning may add a new node here,...
                #assert all(s1 == s2 or self.get_barrier(s1,s2) is not None \
                #        for (s1, s2) in product(basin1, basin2))
                keep |= basin1 | basin2

        for node in list(self.nodes):
            self.nodes[node]['pruned'] = self.nodes.get('pruned', 0)
            if node in keep:
                self.nodes[node]['pruned'] = 0
            else:
                # NOTE: this is potentially problematic!
                # -> now we have two types of inactive nodes ...
                self.nodes[node]['active'] = False
                self.nodes[node]['pruned'] += 1
