# -*- coding: utf-8 -*-
"""DNA-RNA-Transformer landscape exploration.

This module provides the Transformer Landscape Object. There are four basic 
features: expansion, coarse graining, simulation and pruning. 

TrafoLandscape.expand: 
    Identifies parent structures (all sufficiently populated local minima --
    or whatever remains after previous coarse graining, simulation and pruning).
     * The parent structures get connected to the current MFE structure. This
        is a heuristic, that sorts parents by base-pair distance to the MFE, and
        then uses the barrier energy as an upper bound for subsequent searches.
        Eventually, the MFE structure is connected to the parent structure with
        the highest rate into the MFE structure (and potentially to other
        structures with a lower base-pair distance).
    Every parent structure produces and connects its fraying neighbors. 

TrafoLandscape.coarse_grain: 
    Identifies all active structures (parent structures and newly identified
    neighbors after expansion) and merges them using a local flooding
    procedure. The occupancy of a structure is transferred to one (or more)
    energetically better structures, if the barrier toward that structure is
    minimal, and smaller than some given threshold. A merging procedure
    connects all neighbors of the merged node directly to the new local
    minimum.

TrafoLandscape.simulate:
    ...

TrafoLandscape.prune:
    Identifies all insufficiently populated structures (after simulation)
    and merges them using a local flooding procedure. The occupancy of a
    structure is transferred to one (or more) energetically better active
    structures, if the barrier toward that structure is minimal.  A merging
    procedure connects all neighbors of the pruned node directly to the new
    local minimum.

On transition rates:
    * There is always an edge between two nodes unless there has never been the
    attempt to connect two nodes.
    * An edge has saddle energy float('inf') if there exists some path between
    the nodes, but that path has an unknown barrier. The saddle energy can be updated
    if the parameters for findpath are more relaxed than in previous attempts, or if
    an indirect path barrier becomes known.
    * An edge has a valid saddle energy if the energy barrier has been determined
    by a direct or indirect path connecting those structures. A known saddle energy
    cannot be updated with findpath calls.

Coded by: Stefan Badelt <stef@tbi.univie.ac.at>
University of Vienna, Department of Theoretical Chemistry
"""
from __future__ import division, print_function
from builtins import map
from builtins import range

import re
import math
import networkx as nx
import subprocess as sub
from struct import pack
from itertools import combinations, product, islice

import RNA
import ribolands as ril

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Tracks findpath calls for profiling output.                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
PROFILE = {'findpath-calls': 0,
           'fraying1': 0,
           'fraying2': 0,
           'mfe': 0,
           'connect': 0,
           'cogr': 0,
           'prune': 0}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Cache to look-up base-pair distances                                         #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
BPD_CACHE = {}

def get_bpd_cache(s1, s2):
    global BPD_CACHE
    if (s1, s2) in BPD_CACHE:
        return BPD_CACHE[(s1, s2)]
    elif (s2, s1) in BPD_CACHE:
        return BPD_CACHE[(s2, s1)]
    else :
        dist = RNA.bp_distance(s1, s2)
        BPD_CACHE[(s1, s2)] = dist
        return dist

def clear_bpd_cache(TL):
    for k, v in list(BPD_CACHE.items()):
        if not TL.has_edge(*k):
            del BPD_CACHE[k]
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Cache to fold exterior loops - not used at the moment.                       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
EFOLD_CACHE = {}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Custom error definitions                                                     #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
class TrafoUsageError(Exception):
    pass

class TrafoAlgoError(Exception):
    pass

class DebuggingAlert(Exception):
    pass

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Transformer Landscape Object                                                 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
class TrafoLandscape(nx.DiGraph):
    """ 
    A cotranscriptional interface to explore RNA energy landscapes.

    A directed graph (:ob:`networkx.DiGraph()`) where nodes are RNA secondary
    structures and edges are transition rates between secondary structures.  Two
    structures are called neighbors, if there exists an edge connecting them. 

    This object is initialized with the full-length RNA sequence and general RNA
    folding parameters. New (or initial) conformations are found using the
    routine TrafoLandscape.expand(), the landscape can be coarse-grained into a
    more compact representation using the routine TrafoLandscape.coarse_grain()
    and TrafoLandscape.get_simulation_files_tkn() returns a rate matrix to
    simulate RNA folding, TrafoLandscape.update_occupancies_tkn() reads a
    simulation output and updates the occupancies of structures accordingly.
    TrafoLandscape.prune() removes improbable nodes from the landscape.

    Args:
        fullseq (str): The nucleotide sequence of a full-length molecule.
        vrna_md (:obj:`RNA.md()`): ViennaRNA model details. Contains RNA
            folding parameters such as Temperature, noLP, etc...
    """
 
    def __init__(self, fullseq, vrna_md):
        super(TrafoLandscape, self).__init__()

        self._full_sequence = fullseq
        self._model_details = vrna_md
        self._fold_compound = RNA.fold_compound(fullseq, vrna_md)

        # Calculate full backtracking matrix for all subsequent MFE runs.
        _ = self._fold_compound.mfe()

        # Adjust simulation parameters
        self._RT = 0.61632077549999997
        if vrna_md.temperature != 37.0:
            kelvin = 273.15 + vrna_md.temperature
            self._RT = (self._RT/310.15) * kelvin

        # Private instance variables:
        self._transcript_length = 0
        self.total_time = 0
        self._nodeid = 0

        # Default parameters:
        self._fpath = 20    # findpath_search_width

        self._k0 = 2e5    # set directly
        self._dG_max = 0  # set using t_slow
        self._dG_min = 0  # set using t_fast

    @property
    def full_sequence(self):
        return self._full_sequence

    @property
    def transcript(self):
        return self._full_sequence[0:self._transcript_length]

    @property
    def t_fast(self):
        """Time-scale separation parameter for coarse-graining.

        Expected folding times faster than t-fast are considered instantaneous,
        and used to assign a conformation to a macrostate. This parameter is
        equivalent to the internal dG_min parameter, which quantifies the minimal
        energy barrier to separate two macrostates.

        t_fast = 1/(self._k0 * exp(-self.dG_min/self._RT))
        """
        return 1/(self._k0 * math.exp(-self._dG_min/self._RT))

    @t_fast.setter
    def t_fast(self, value):
        self._dG_min = max(0, -self._RT * math.log(1 / value / self._k0))

    @property
    def t_slow(self):
        """Time-scale separation parameter for finding new neighbors.

        When new a potential new relevant RNA conformation was found, then this
        conformation needs to be connected to the existing energy landscape.
        This parameter is used to *reject* a transition rate toward a new
        energetically better structure if it is lower than t-slow, it is
        equivalent to the internal dG_max parameter, which quantifies the maximal
        energy barrier that can be overcome within a folding simulation.

        t_slow = 1(self._k0 * exp(-self.dG_max/self._RT))
        """
        if self._dG_max == 0:
            return None
        else:
            return 1/(self._k0 * exp(-self._dG_max/self._RT))

    @t_slow.setter
    def t_slow(self, value):
        if value is None:
            self._dG_max = 0
        else:
            self._dG_max = -self._RT * math.log(1 / value / self._k0)

    @property
    def findpath_search_width(self):
        """Search width for a heuristic to find transition state energies.
        """
        return self._fpath

    @findpath_search_width.setter
    def findpath_search_width(self, val):
        self._fpath = val

    def graph_copy(self):
        """Returns a copy of the TrafoLandscape Graph. 

        This does not include internal parameters such as the current
        transcript length, current node ID, current time, etc.
        """
        copy = TrafoLandscape(self.full_sequence, self._model_details)
        copy.add_nodes_from(self.nodes(data = True))
        copy.add_edges_from(self.edges(data = True))
        return copy

    def graph_to_json(self, name):
        """Prints the current graph into a JSON file format.

        There exits a script to visualize the output files using d3js in your
        browser. Search for it in: ribolands/d3js/start_server.py

        Args: 
            name (str): The basename of the *.json file.
        """
        import json
        from networkx.readwrite import json_graph
        d = json_graph.node_link_data(self)
        json.dump(d, open(name+'.json', 'w'))

    def graph_to_pdf(self, name):
        import matplotlib.pyplot as plt
        name += 'mpl.pdf'
        nl = [x for x in self.nodes if self.nodes[x]['active']]
        nd = nx.get_node_attributes(self, 'identity')
        nx.drawing.nx_pylab.draw_networkx(self, nodelist=nl, labels=nd)

        plt.axis('off')  # turn of axis
        plt.savefig(name)
        plt.clf()

    def sorted_nodes(self, descending = False):
        """ Returns active nodes and their attributes sorted by energy.

        Args:
          descending (bool, optional): sorting parameter.
            True: energetically high to low.
            False: energetically low to high.
            Defaults to False.

        """
        active = [n_d for n_d in self.nodes(data = True) if n_d[1]['active']]
        return sorted(active, key = lambda x: (
            x[1]['energy'], x[0]), reverse = descending)

    def has_active_edge(self, s1, s2):
        return self.has_edge(s1, s2) and self.get_saddle(s1, s2) < float('inf')

    def get_saddle(self, s1, s2):
        """Returns the saddle energy of a transition edge."""
        if self.has_edge(s1, s2):
            return self[s1][s2]['saddle']
        else:
            return None

    def get_barrier(self, s1, s2):
        """Returns the barrier energy of a transition edge."""
        if self.has_edge(s1, s2):
            return self[s1][s2]['saddle'] - self.nodes[s1]['energy']
        else:
            return None

    def get_rate(self, s1, s2):
        """Returns the direct transition rate of two secondary structures."""
        if self.has_edge(s1, s2):
            return self[s1][s2]['weight']
        else:
            return 0

    def add_transition_edges(self, s1, s2, ts = None, 
            maxbar = float('inf'), fpathW = None, call = None):
        """Calculates transition rates from (multiple) direct path barrier heights.

        Uses the *findpath* direct path heuristic to find the lowest energy
        barrier between all local features of two secondary structures.
        Typically s2 is the new, energetically better structure, but this is
        not enforced. In case a transient structure "ts" is specified, the rate
        is calculated as the minimum between the direct path s1 -> s2 and the
        indirect folding path s1 -> ts -> s2.  Rates are computed using the
        Arrhenius model with global parameters k0 and RT.

        If this function finds a valid saddle energy between s1 and s2, then it adds 
        all necessary nodes along the folding path to the graph.
        
        Args:
          s1 (str): start secondary structure (must be part of TrafoLandscape already)
          s2 (str): final secondary structure (may be added to TrafoLandscape)
          ts (str, optional): transient secondary structure which is not on the
            direct folding path. If a transient structure is specified, the direct
            path barriers s1 -> ts and ts -> s2 must be known already.
          maxbar (flt, optional): an upper bound on the energy barrier from s1 -> s2.
            Defaults to float('inf'). 
          fpathW (int, optional): the search width parameter for the findpath routine.
            Defaults to None: 2 * base-pair distance.

        Returns:
          bool: True if a non-zero rate between structures has been found.
        """
        assert s1 != s2
        assert self.has_node(s1)
        if maxbar < self._dG_min:
            maxbar = self._dG_min

        def findpath_wrap(s1, s2, maxbar = float('inf'), fpathW = None):
            # Returns the maximum barrier for the overall refolding.
            features = get_bp_change(fullseq, s1, s2)
            s1en = self.nodes[s1]['energy']

            totbar = 0
            for [cseq, cs1, cs2] in features:
                fpw = 2 * get_bpd_cache(cs1, cs2) if fpathW is None else fpathW
                barrier, gain, _, _ = get_fpath_cache(cseq, cs1, cs2, md, maxbar, fpw)
                totbar += barrier
                #if barrier > limiting_barrier:
                #    limiting_barrier = barrier

            return round(self.nodes[s1]['energy'] + totbar, 2)

        fullseq = self._full_sequence
        md = self._model_details
        fc = self._fold_compound

        _RT = self._RT
        _k0 = self._k0

        if ts: # Lookup the in-direct path barrier first
            tsE1 = self.get_saddle(s1, ts)
            tsE2 = self.get_saddle(ts, s2)
            tsE = max(tsE1, tsE2)
            assert tsE != float('inf')
            # Just to be sure, because right now there should be no situation
            # where maxbar and ts are both non-default. You may very well
            # remove this statement if that changes.
            assert maxbar == float('inf')
        elif self._dG_max: 
            maxbar = min(maxbar, self._dG_max)  

        saddleE = self.get_saddle(s1, s2)
        #print('se', saddleE)
        if saddleE is None or saddleE == float('inf'):
            saddleE = findpath_wrap(s1, s2, maxbar, fpathW)
            #print('se2', saddleE)
            # Collect some profiling data.
            PROFILE['findpath-calls'] += 1
            if call: PROFILE[call] += 1
        if ts and tsE < saddleE:
            #print('ts', tsE)
            saddleE = tsE

        assert saddleE is not None

        if saddleE == float('inf') and self.has_node(s2):  # Add the edge.
            self.add_weighted_edges_from([(s1, s2, 0)])
            self.add_weighted_edges_from([(s2, s1, 0)])
            self[s1][s2]['saddle'] = float('inf')
            self[s2][s1]['saddle'] = float('inf')
        elif saddleE < float('inf'):  # Add the edge.
            e1 = self.nodes[s1]['energy']
            e2 = self.nodes[s2]['energy'] if self.has_node(s2) \
                                         else round(fc.eval_structure(s2), 2)

            # ensure saddle is not lower than s1, s2
            saddleE = max(saddleE, max(e1, e2))

            # Energy barrier
            dG_1s = saddleE - e1
            dG_2s = saddleE - e2

            #if s1 == '.......(((((((...)))))))((((((.......))).)))......':
            #    print(s1, e1, saddleE, dG_1s)
            #    print(s2)
            #elif s2 == '.......(((((((...)))))))((((((.......))).)))......':
            #    print(s1, e1, saddleE)
            #    print(s2, e2, dG_2s)
            # Arrhenius model
            k_12 = _k0 * math.exp(-dG_1s/_RT)
            k_21 = _k0 * math.exp(-dG_2s/_RT)
            self.add_weighted_edges_from([(s1, s2, k_12)])
            self.add_weighted_edges_from([(s2, s1, k_21)])
            self[s1][s2]['saddle'] = saddleE
            self[s2][s1]['saddle'] = saddleE

        return saddleE != float('inf')

    def expand(self, extend = 1, exp_mode = 'default', mfree = 6):
        """Find new secondary structures and add them to :obj:`TrafoLandscape()`

        The function supports two move-sets: 1) The mfe structure for the current
        sequence length is connected to all present structures, 2) The conformation
        graph is expanded using helix-fraying.

        Args:
          extend (int, optional): number of nucleotide extensions before graph
            expansion (updates the global variable transcript length). Defaults to 1.
          exp_mode (str, optional): choose from "mfe-only": only use current mfe
            structure as potential new neighbor. "fraying-only": only use fraying 
            neighborhood. "default": do both mfe and fraying.
          mfree (int, optional): minimum number of freed bases during a
            helix-opening step. Defaults to 6.

        Returns:
          int: Number of new nodes
        """
        fseq = self.full_sequence
        self._transcript_length += extend
        if self._transcript_length > len(fseq):
            self._transcript_length = len(fseq)
        seq = self.transcript

        csid = self._nodeid
        md = self._model_details
        fc = self._fold_compound

        if exp_mode not in ('default', 'fullconnect'):
            raise TrafoUsageError('Invalid expansion mode used!')

        # Calculate MFE of current transcript
        mfess, mfe = fc.backtrack(len(seq))
        future = '.' * (len(fseq) - len(seq))
        mfess = mfess + future

        # If there is no node because we are in the beginning, add the node.
        if len(self) == 0:
            en = round(fc.eval_structure(mfess), 2)
            self.add_node(mfess, energy = en, occupancy = 1.0, 
                          identity = self._nodeid, active = True, last_seen = 0)
            self._nodeid += 1

        # Save the set of active parent nodes for later.
        parents = [n for n, d in self.nodes(data = True) if d['active']]

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # Connect MFE to the parent ensemble.                                  #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        num = self.expand_connect_mfe(parents, mfess, ddG = float('inf'), exp_mode = exp_mode)
        if num == 0 and len(parents) > 1:
            assert (not self.has_node(mfess) or \
                (self.has_node(mfess) and not self.nodes[mfess]['active']))
            print("# WARNING: mfe secondary structure not connected\n# {}".format(
                    mfess[0:self._transcript_length]))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # Fraying neighbors (1/2): find and connect new structures to parents. #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

        new_nodes = self.expand_fraying_neighbors(parents, mfree = 6)

        # Post processing of graph after expansion:
        # remove nodes that have been inactive for a long time.
        for ni in list(self.nodes()):
            if self.nodes[ni]['active']:
                self.nodes[ni]['last_seen'] = 0
            else:
                self.nodes[ni]['last_seen'] += 1
            if self.nodes[ni]['last_seen'] >= 5:
                self.remove_node(ni)

        if exp_mode == 'fullconnect':
            # Just connect everything to everything and be done with it.
            nlist = self.sorted_nodes()
            for ((ni, di), (nj, dj)) in combinations(nlist, 2):
                _ = self.add_transition_edges(ni, nj, call='connect')
                assert _ # assert statements may be disabled during program execution.
            return self._nodeid - csid

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # Fraying neighbors (2/2): connect neighbors of neighboring parents.   #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

        for p1, p2 in combinations(new_nodes.keys(), 2):
            if not self.has_active_edge(p1, p2): continue
            for (np1, np2) in product(new_nodes[p1], new_nodes[p2]):
                if np1 == np2: continue
                dd = get_bpd_cache(np1, np2)
                d1 = get_bpd_cache(np1, p1)
                d2 = get_bpd_cache(p1, p2)
                d3 = get_bpd_cache(p2, np2)
                if dd <= max(d1, d2, d3):
                    #TODO: check that, changed to barrier
                    tsE1 = self.get_barrier(np1, p1)
                    tsE2 = self.get_barrier(p1, p2)
                    tsE3 = self.get_barrier(p2, np2)
                    tsE = max(tsE1, tsE2, tsE3)
                    assert tsE != float('inf')
                    _ = self.add_transition_edges(np1, np2, 
                            maxbar = round(tsE + 0.10, 2), call = 'fraying2')

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # Triangle connect:                                                    #
        # -----------------                                                    #
        # For each node, if two neighbors are not connected, and if the        #
        # basepair distance between neighbors is shorter than one of the       #
        # distances to the neighbor, then try connect them with upperbound.    #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
        update = True
        while (update):
            update = False
            for np in parents:
                nbrs = [n for n in self.successors(np) if self.nodes[n]['active'] and \
                        self.has_active_edge(n, np)]

                for n1, n2 in combinations(nbrs, 2):
                    #if self.has_active_edge(n1, n2):
                    if self.has_edge(n1, n2): 
                        # If it is inactive, that means that a better, indirect
                        # connection has been found.
                        continue
                    dd = get_bpd_cache(n1, n2)
                    d1 = get_bpd_cache(np, n1)
                    d2 = get_bpd_cache(np, n2)
                    if dd <= max(d1, d2):
                        tsE1 = self.get_barrier(np, n1)
                        tsE2 = self.get_barrier(np, n2)
                        tsE = max(tsE1, tsE2)
                        assert tsE != float('inf')
                        _ = self.add_transition_edges(n1, n2, 
                                maxbar = round(tsE + 0.10, 2), call = 'connect')
                        update = update or _

        clear_bpd_cache(self)

        return self._nodeid - csid

    def expand_connect_mfe(self, parents, mfess, ddG = float('inf'), exp_mode = 'default'):
        """ Connect the MFE structure to parent structures.

        This function checks the node attibutes:
            TL.node.data: active, energy
        and sets:
            TL.node.data: active, last_seen
    
        Args:
            parents (list): List of secondary structures, they must be 'active'.
            mfess (str): The MFE structure.
            ddG (float): An upper bound on the barrier energy from parents to MFE
    
        Returns: 
          int: The number of active connections from/to the MFE structure.
        """
        TL = self
        fc = TL._fold_compound
    
        connections = 0
        for parent in sorted(parents, key = lambda x: get_bpd_cache(mfess, x)):
            if parent == mfess: continue
            assert TL.nodes[parent]['active']
    
            # If there is an active edge, make sure the MFE becomes an active node.
            if TL.has_active_edge(parent, mfess):
                TL.nodes[mfess]['active'] = True
                TL.nodes[mfess]['last_seen'] = 0
    
            # If there is a node, but no active edge, try to add it and activate the node.
            elif TL.has_node(mfess):
                if TL.add_transition_edges(parent, mfess, 
                        maxbar = ddG + 0.01, call = 'mfe'):
                    # in case it was there but inactive
                    TL.nodes[mfess]['active'] = True
                    TL.nodes[mfess]['last_seen'] = 0
    
            # If there is no node, try to add a new edge & initialize the node.
            elif TL.add_transition_edges(parent, mfess, 
                        maxbar = ddG + 0.01, call = 'mfe'):
                    en = round(fc.eval_structure(mfess), 2)
                    TL.nodes[mfess]['active'] = True
                    TL.nodes[mfess]['last_seen'] = 0
                    TL.nodes[mfess]['energy'] = en
                    TL.nodes[mfess]['occupancy'] = 0.0
                    TL.nodes[mfess]['identity'] = TL._nodeid
                    TL._nodeid += 1
    
            if exp_mode == 'default' and TL.has_active_edge(mfess, parent):
                barrier = TL.get_saddle(mfess, parent) - TL.nodes[parent]['energy']
                barrier = max(barrier, TL._dG_min)
                if barrier < ddG:
                    ddG = barrier
                connections += 1
            elif exp_mode == 'fullconnect':
                assert TL.has_active_edge(mfess, parent)
                connections += 1
    
        return connections

    def expand_fraying_neighbors(self, parents, mfree = 6):
        """Eyoo -- test me!"""
        seq = self.transcript
        fc = self._fold_compound
        md = self._model_details

        # This dictionary is a cache for feature expansion:
        #   ext_moves[ext_seq] = [structure, set((con, paren), ...)]
        # where ext_seq = exterior-loop sequence with ((xxx)) constraints
        ext_moves = dict()

        # Keep track of all new nodes to connect them with each other.
        new_nodes = dict() 

        def ext_move_speedup(ext_seq):
            """ Connect the neighbor with *historic* transitions of parents.

            We store the exterior-open neighbor here, that means there are three
            possible reasons for duplication:
              1) different (or longer) helix was opened / same historic features
              2) the same helix was opened / difference is in historic features
              3) different helix / different history
            """
            for (parent, child) in ext_moves[ext_seq][1]:
                assert parent != ni  # Parents may never be the same
                if child == nbr:
                    # the parents differ in fraying helices, 
                    # no historic differences
                    continue

                if not self.has_active_edge(parent, ni):
                    # no need to worry about it.
                    continue

                # Now the case with historic differences ...
                assert self.has_node(child) and self.has_node(nbr)
                assert self.nodes[child]['active'] and self.nodes[nbr]['active']

                if self.has_active_edge(nbr, child):
                    continue

                # Calculate saddleE from saddleE of parents!
                # The forward barrier parent -> saddle and reverse barrier ni -> saddle
                # are added to the new structures child -> sp and nbr -> sp
                sp = self.get_saddle(parent, ni)
                fwbar = sp - self.nodes[parent]['energy'] 
                bwbar = sp - self.nodes[ni]['energy'] 
                fwsp = round(self.nodes[child]['energy'] + fwbar, 2)
                bwsp = round(self.nodes[nbr]['energy'] + bwbar, 2)

                if fwsp == bwsp:  # avoid findpath!
                    self.add_weighted_edges_from([(child, nbr, None)])
                    self.add_weighted_edges_from([(nbr, child, None)])
                    self[child][nbr]['saddle'] = fwsp
                    self[nbr][child]['saddle'] = bwsp
                    if self.add_transition_edges(nbr, child):
                        pass
                    else:
                        raise TrafoAlgoError('Did not add historic edge!')
            return

        for ni in parents:
            new_nodes[ni] = []

            # short secondary structure (without its future)
            sss = ni[0:len(seq)]

            # compute a set of all helix fraying open steps
            opened = open_fraying_helices(seq, sss, mfree)

            # Do a constrained exterior loop folding for all fraying structures
            # and connect them to the parent conformation.
            connected = set([ni]) # add the parent
            for onbr in opened:
                nbr, ext_seq = fold_exterior_loop(md, seq, onbr, ext_moves)

                future = '.' * (len(ni) - len(nbr))
                nbr += future

                # Not a new structure
                if nbr in connected: 
                    continue
                connected.add(nbr)

                # Add fraying helix transitions 
                if self.has_active_edge(ni, nbr):
                    # NOTE: this may activate a previously coarse-grained node ...
                    # That's not pretty, but at least it is not wrong.
                    self.nodes[nbr]['active'] = True
                    self.nodes[nbr]['last_seen'] = 0
                elif self.has_node(nbr):
                    if self.add_transition_edges(ni, nbr, call = 'fraying1'):
                        self.nodes[nbr]['active'] = True
                        self.nodes[nbr]['last_seen'] = 0
                elif self.add_transition_edges(ni, nbr, call = 'fraying1'):
                    enbr = round(fc.eval_structure(nbr), 2)
                    self.nodes[nbr]['energy'] = enbr
                    self.nodes[nbr]['active'] = True
                    self.nodes[nbr]['last_seen'] = 0
                    self.nodes[nbr]['occupancy'] = 0.0
                    self.nodes[nbr]['identity'] = self._nodeid
                    self._nodeid += 1

                if self.has_node(nbr) and self.nodes[nbr]['active']:
                    #TODO: alright this is fucked up, new nodes can be parents... shit
                    #if nbr not in parents:
                    new_nodes[ni].append(nbr)

                    # This is a shortcut to pre-calculate an edge from the next step:
                    # Every identical ext-change will be connected, if the parents
                    # were connected.
                    if ext_moves[ext_seq][1]:
                        ext_move_speedup(ext_seq)
                    ext_moves[ext_seq][1].add((ni, nbr))
                    assert self.get_barrier(ni,nbr) is not None
                    assert self.get_barrier(ni,nbr) != float('inf')

        return new_nodes

    def coarse_grain(self, dG_min = None):
        """Landscape coarse-graining base on energy barriers.

        Every structure gets assigned to an energetically better or equal
        macrostate. If a structure gets a assigned to itself, well then it
        represents the macrostate

        Try it in two iterations:
            -> assign every structure to its macrostate. 
            -> remove nodes but and use transition to 

        Processes an energetically sorted list of structures (high to low
        energies) and tries to merge their occupancy into a neighboring, better
        conformation (forming a macro-state). The current structure is merged
        (and therefore removed from the graph) if there exists a neighbor that
        has better energy and the energy barrier is lower than the dG_min
        parameter. If there are multiple better neighbors with a minimal
        transition energy barrier, then the occupancy is transfered to the
        neighbor with the lowest energy, in the degenerate case the
        lexicographically first structure is chosen.

        Args:
          dG_min (flt, optional): Minimum energy barrier between separate
            macrostates.

        Returns:
          dict[del-node] = macro-node: A mapping from deleted nodes to macro-state
        """

        merged_nodes = dict()
        merged_to = dict()

        if dG_min is None:
            dG_min = self._dG_min + 0.01
        assert dG_min is not None

        # sort by energy (high to low)
        for ni, data in sorted(self.nodes(data = True),
                               key = lambda x: (x[1]['energy'], x), reverse = True):
            #print(ni, data)

            if data['active'] == False:
                continue
            en = data['energy']

            # get all active neighbors (low to high) including the high
            # neighbors, bec we have to connect them later
            nbrs = [x for x in sorted(self.successors(ni), 
                key = lambda y: (self.nodes[y]['energy'], y), 
                    reverse = False) if self.nodes[x]['active']]

            nbrs = [x for x in nbrs if self.get_saddle(ni, x) != float('inf')]

            if nbrs == []:
                break

            # lowest neighbor structure and energy
            best, been = nbrs[0], self.nodes[nbrs[0]]['energy']

            if been - en > 0.0001:
                # local minimum
                continue

            # among all energetically equal and lower neighbors, find the
            # neighbor(s) with the lowest energy barrier ...
            (transfer, minsE) = (set([best]), self.get_saddle(ni, best))
            for e, nbr in enumerate(nbrs[1:]):
                if self.nodes[nbr]['energy'] - en >= 0.0001:
                    break
                sE = self.get_saddle(ni, nbr)
                if sE + 0.01 < minsE : # a truly lower energy saddle
                    (transfer, minsE) = (set([nbr]), sE)
                elif abs(sE - minsE) < 0.0001 : # an equal conformation.
                    transfer.add(nbr)
                    minsE = sE

            if minsE - en - dG_min > 0.0001:  # avoid precision errors
                # do not merge, if the barrier is too high.
                continue

            # connect all neighboring nodes to the transfer nodes
            for nb1 in nbrs:
                for trans in transfer:
                    if nb1 == trans:
                        continue
                    if self.nodes[nb1]['energy'] > self.nodes[trans]['energy']:
                        (s1, s2) = (nb1, trans) 
                    else:
                        (s1, s2) = (trans, nb1) 
                    always_true = self.add_transition_edges(s1, s2, ts = ni, call = 'cogr')
                    if always_true is False:
                        print(s1)
                        print(ni)
                        print(s2)
                        raise TrafoAlgoError('Did not add the transition edge!')

            # remove the node
            self.nodes[ni]['active'] = False
            self.nodes[ni]['last_seen'] = 1
            for trans in transfer:
                self.nodes[trans]['occupancy'] += self.nodes[ni]['occupancy']/len(transfer)
            self.nodes[ni]['occupancy'] = 0.0

            # TODO: need to double check this...
            merged_nodes[ni] = transfer
            for trans in transfer:
                if trans in merged_to:
                    merged_to[trans].append(ni)
                else:
                    merged_to[trans] = [ni]

            if ni in merged_to:
                fathers = merged_to[ni]
                for f in fathers:
                    merged_nodes[f] = transfer
                    for trans in transfer:
                        merged_to[trans].append(f)
                del merged_to[ni]

        #for (k,v) in merged_nodes.items():
        #    print(k, v)
        return merged_nodes

    def simulate(self, t0, t8, tmpfile = None):
        # treekin wrapper function using:
        #   "self.get_simulation_files_tkn"
        #   "self.update_occupancies_tkn"
        raise NotImplementedError

    def get_simulation_files_tkn(self, name, binrates = True):
        """ Print a rate matrix and the initial occupancy vector.

        This function prints files and parameters to simulate dynamics using the
        commandline tool treekin. A *.bar file contains a sorted list of present
        structures, their energy and their neighborhood and the corresponding
        energy barriers. A *.rts or *.rts.bin file contains the matrix of
        transition rates either in text or binary format. Additionaly, it returns
        a vector "p0", which contains the present occupancy of structures. The
        order or elements in p0 contains

        Note:
          A *.bar file contains the energy barriers to transition between local
          minima. In contrast to files produced by `barriers`, where local minimum
          is always *directly* connected to an energetically better local minimum,
          here a path towards the MFE structure can proceed via an energetically
          worse structure first.

        Args:
          name (str): Name of output files name.bar, name.rts, name.rts.bin.
          binrates (bool, optional): Print rates in binary format or text format.
            Defaults to True: binary format.

        """
        seq = self.transcript

        sorted_nodes = self.sorted_nodes(descending = False)
        num = len(sorted_nodes)+1

        bfile = name + '.bar'
        rfile = name + '.rts'
        brfile = rfile + '.bin'
        p0 = []

        with open(bfile, 'w') as bar, open(rfile, 'w') as rts, open(brfile, 'wb') as brts:
            bar.write("  ID {}  Energy  {}\n".format(seq, 
                ' '.join(map("{:7d}".format, range(1,num)))))
            brts.write(pack("i", len(sorted_nodes)))
            for e, (ni, data) in enumerate(sorted_nodes, 1):
                # Calculate barrier heights to all other basins.
                mystr = ''
                for ee, (be, _) in enumerate(sorted_nodes, 1):
                    if e == ee:
                        #mystr += '       '
                        mystr += ' {:7.2f}'.format(float('nan'))
                        continue
                    sE = self.get_saddle(be, ni)
                    if sE is not None:
                        mystr += ' {:7.2f}'.format(sE - data['energy'])
                    else:
                        #mystr += '       '
                        mystr += ' {:7.2f}'.format(float('nan'))

                # Print structures and neighbors to bfile:
                bar.write("{:4d} {} {:7.2f} {}\n".format(
                    e, ni[:len(seq)], data['energy'], mystr))

                # Add ni occupancy to p0
                if data['occupancy'] > 0:
                    p0.append("{}={}".format(e, data['occupancy']))

                # Print rate matrix to rfile and brfile
                trates = []
                rates = []
                for (nj, jdata) in sorted_nodes:
                    if self.has_edge(ni, nj):
                        rates.append(self[ni][nj]['weight'])
                        trates.append(self[nj][ni]['weight'])
                    else:
                        rates.append(0)
                        trates.append(0)
                line = "".join(map("{:10.4g}".format, rates))
                rts.write("{}\n".format(line))
                for r in trates:
                    brts.write(pack("d", r))

        return [bfile, brfile if binrates else rfile, p0, sorted_nodes]

    def load_barfile(self, bfile):
        """Write a function that initializeds the graph using a present
        bar-file. We can use that to speed up the second time where we simulate
        just for different visualization (e.g.  plot-minh, svg, etc.)
        """
        pass

    #   update_time_and_occupancies_tkn(self, tfile)
    def update_occupancies_tkn(self, tfile, sorted_nodes):
        """
          Update the occupancy in the Graph and the total simulation time
        """
        # http://www.regular-expressions.info/floatingpoint.html
        reg_flt = re.compile(b'[-+]?[0-9]*.?[0-9]+([eE][-+]?[0-9]+)?.')

        lastlines = sub.check_output(['tail', '-2', tfile]).strip().split(b'\n')
        if not reg_flt.match(lastlines[0]):
            raise TrafoAlgoError('Cannot parse simulation output', tfile)
        else:
            if reg_flt.match(lastlines[1]):
                time = float(lastlines[1].split()[0])
                iterations = None
                tot_occ = sum(map(float, lastlines[1].split()[1:]))
                for e, occu in enumerate(lastlines[1].split()[1:]):
                    ss = sorted_nodes[e][0]
                    self.nodes[ss]['occupancy'] = float(occu)/tot_occ
            else :
                time = float(lastlines[0].split()[0])
                iterations = int(lastlines[-1].split()[-1])
                tot_occ = sum(map(float, lastlines[0].split()[1:]))
                for e, occu in enumerate(lastlines[0].split()[1:]):
                    ss = sorted_nodes[e][0]
                    self.nodes[ss]['occupancy'] = float(occu)/tot_occ

        return time, iterations

    def prune(self, p_min, detailed = True, keep_reachables = False):
        """ Delete nodes or report them as still reachable.

        Use the occupancy cutoff to choose which nodes to keep and which ones to
        remove. Every node with occupancy < cutoff will be removed and its neighbors
        connected with each other. 

        Arguments:
          p_min (flt): Occupancy cutoff for neighbor generation.
            Defaults to None: using global TrafoLandscape parameter.

        Returns:
          int, int, int:
            number of deleted nodes,
            number of still reachable nodes,

        """
        deleted_nodes = 0
        still_reachables = 0

        for ni, data in self.sorted_nodes(descending = False):  # sort high to low..
            if data['occupancy'] >= p_min:
                continue
            en = data['energy']

            # get all active neighbors (low to high)
            nbrs = [x for x in sorted(self.successors(ni), 
                key=lambda x: self.nodes[x]['energy'], reverse = False) \
                        if self.nodes[x]['active']]

            nbrs = [x for x in nbrs if self.get_saddle(ni, x) != float('inf')]

            if nbrs == []:
                raise TrafoAlgoError('there must be some valid neighbors! :-/')

            # lowest neighbor structure and energy
            best, been = nbrs[0], self.nodes[nbrs[0]]['energy']

            if keep_reachables and been - en > 0.0001:
                raise DeprecationWarning('do not keep keep reachables!')
                still_reachables += 1
                continue

            # among *all* neighbors, find the neighbor with the lowest energy barrier
            (transfer, minsE) = (best, self.get_saddle(ni, best))
            for e, nbr in enumerate(nbrs[1:]):
                sE = self.get_saddle(ni, nbr)
                if sE - minsE < 0.0001:
                    (transfer, minsE) = (nbr, sE)

            # remove the node
            self.nodes[ni]['active'] = False
            self.nodes[ni]['last_seen'] = 1
            self.nodes[transfer]['occupancy'] += self.nodes[ni]['occupancy']
            self.nodes[ni]['occupancy'] = 0.0
            deleted_nodes += 1

            for e, nb1 in enumerate(nbrs, 1):
                for nb2 in nbrs[e:]:
                    assert detailed
                    #NOTE: we could switch between detailed and fake:
                    #   - if we are connecting to the bp-dist closest one: detailed
                    #   - if we are connecting to sthg with shorter bp-dist: detailed
                    always_true = self.add_transition_edges(nb2, nb1, ts = ni, call = 'prune')
                    if always_true is False:
                        raise TrafoAlgoError('Did not add the transition edge!')

        return deleted_nodes, still_reachables

    def sorted_trajectories_iter(self, sorted_nodes, tfile, softmap = None):
        """ Yields the time course information using a treekin output file.

        Args:
          sorted_nodes (list): a list of nodes sorted by their energy
          tfile (str): treekin-output file name.
          softmap (dict, optional): A mapping to transfer occupancy between
            states. Likely not the most efficient implementation.

        Yields:
          list: ID, time, occupancy, structure, energy
        """
        # http://www.regular-expressions.info/floatingpoint.html
        reg_flt = re.compile('[-+]?[0-9]*.?[0-9]+([eE][-+]?[0-9]+)?.')

        ttime = self.total_time

        with open(tfile) as tkn:
            # this is ugly, but used to check if we're at the last line
            prevcourse = []
            tknlines = tkn.readlines()
            for line in tknlines:
                if reg_flt.match(line):
                    course = list(map(float, line.strip().split()))
                    time = course[0]

                    # softmap hack:
                    # preprocess the timeline by merging all states
                    if softmap:
                        macrostates = [0] * len(course)
                        macromap = dict()
                        for e, occu in enumerate(course[1:]):
                            ss = sorted_nodes[e][0]

                            # map occupancy to (energetically better)
                            if ss in softmap:
                                mapss = softmap[ss]     # set of mapss
                                mapids = [macromap[ma] for ma in mapss]
                            else:
                                # we *must* have seen this state before, given
                                # there are no degenerate sorting errors...
                                mapids = [e + 1]
                                macromap[ss] = mapids[0]

                            for mapid in mapids:
                                macrostates[mapid] += occu/len(mapids)

                        course[1:] = macrostates[1:]

                    for e, occu in enumerate(course[1:]):
                        # is it above visibility threshold?
                        ss = sorted_nodes[e][0]
                        sss = ss[0:self._transcript_length]

                        yield self.nodes[ss]['identity'], ttime + time, occu, \
                            sss, self.nodes[ss]['energy']
                    prevcourse = course
        return


def open_fraying_helices(seq, ss, free = 6):
    """Generate structures with opened fraying helices. 
    
    Fraying helices share a base-pair with an exterior loop region. The
    function returns n structures, where the first n-1 correspond to individual
    fraying helices opened respectively and structure n has all fraying helices
    opened at once.

    Args:
        seq (str): Primary structure.
        ss (str): Secondary structure.
        free (int, optional): Minimal number of bases freed by a fraying helix
            move. If less bases are freed, and there exists a nested helix, then
            that helix is opened as well. Defaults to 6.
    """
    nbrs = set()
    pt = ril.make_pair_table(ss, base = 0)

    # mutable secondary structure
    nbr = list(ss)

    rec_fill_nbrs(nbrs, ss, nbr, pt, (0, len(ss)), free)

    nbrs.add(''.join(nbr))

    return nbrs

def rec_fill_nbrs(nbrs, ss, mb, pt, myrange, free):
    """Recursive helix opening

    TODO: Function needs testing, but looks good.

    Args:
        nbrs (set): A set of all neighboring conformations.
        ss (str): Reference secondary structure
        mb (list): A mutable version of ss, which will have all fraying helices
            opened at once.
        pt (list): A pair table (zero based)
        myrange (int, int): the range (n, m) of the pt under current
            investigation.  
        free: number of bases that should be freed

    Returns:
        None: mutates the mutable arguments.
    """
    (n, m) = myrange
    skip = 0  # fast forward in case we have deleted stuff
    for i in range(n, m):
        j = pt[i]
        if j == -1:
            continue
        if i < skip:
            continue

        nb = list(ss)
        [o, l] = [0, 0]
        [p, q] = [i, j]

        add = True
        while p < q and (l == 0 or o < free):
            if pt[p] != q or p != pt[q]:
                """ this is a multiloop """
                rec_fill_nbrs(nbrs, ''.join(nb), mb, pt, (p, q), free - o)
                add = False
                break

            # remove the base-pairs
            pt[p] = pt[q] = -1
            nb[p] = nb[q] = '.'
            mb[p] = mb[q] = '.'
            o += 2  # one base-pair deleted, two bases freed

            l = 0  # reset interior-loop size
            while (p < q and pt[p + 1] == -1):
                [p, l] = [p + 1, l + 1]
            p += 1
            while (p < q and pt[q - 1] == -1):
                [q, l] = [q - 1, l + 1]
            q -= 1
            o += l

        if add:
            nbrs.add(''.join(nb))
        skip = j + 1

    return

def fold_exterior_loop(md, seq, con, ext_moves, spacer = 'NNN'):
    """ Constrained folding of the exterior loop.

    All constrained helices are replaced with the motif:
        NNNNNNN
        ((xxx))
    for example a helix with the closing-stack CG-UG:
        CG ~ UG -> CGNNNUG
        (( ~ )) -> ((xxx))
    This reduces the sequence length (n) and therefore the runtime O(n^3),
    and it enables the identification of independent structures with the same
    exterior loop features.

    Args:
        md (RNA.md()):      ViennaRNA model details (temperature, noLP, etc.)
        seq (str):          RNA sequence
        con (str):          RNA structure constraint
        ext_moves (dict()): Dictionary storing all mappings from exterior-loop
                            constraints (features) to parents.
        spacer (str, optional): Specify the sequence of the spacer region.
                            Defaults to 'NNN'.

    Returns:
      (str, str):
    """
    pt = ril.make_pair_table(con, base = 0)
    ext_seq = ''
    ext_con = ''

    # shrink the sequences
    skip = 0
    for i, j in enumerate(pt):
        if i < skip:
            continue
        if j == -1:
            ext_seq += seq[i]
            ext_con += '.'
        else:
            ext_seq += seq[i] + seq[i + 1]
            ext_seq += spacer
            ext_seq += seq[j - 1] + seq[j]
            ext_con += '(('
            ext_con += 'x' * len(spacer)
            ext_con += '))'
            skip = j + 1

    # If we have seen this exterior loop before, then we don't need to
    # calculate again, and we have to trace back if the parents are connected.
    if ext_seq in ext_moves:
        css = ext_moves[ext_seq][0]
    else:
        fc_tmp = RNA.fold_compound(ext_seq, md)
        fc_tmp.constraints_add(
            ext_con, RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP)
        css, cfe = fc_tmp.mfe()
        ext_moves[ext_seq] = [css, set()]
        del fc_tmp

    # replace characters in constraint
    c, skip = 0, 0
    for i, j in enumerate(pt):
        if i < skip:
            continue
        if j == -1:
            con = con[:i] + css[c] + con[i + 1:]
            c += 1
        else:
            c += len(spacer) + 4
            skip = j + 1
    ss = con

    return ss, ext_seq

