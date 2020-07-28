#
# ribolands.trafo
# 
# Home of the TrafoLandscape object.
#

import logging
rlog = logging.getLogger(__name__)

import re
import math
import networkx as nx
import subprocess as sub
from struct import pack
from itertools import combinations, product, islice

import RNA
from ribolands import PrimePathLandscape
from ribolands.utils import make_pair_table
from ribolands.pathfinder import (show_flooded_prime_path, 
                                  BPD_CACHE, get_bpd_cache)

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

def clear_bpd_cache(TL):
    for k, v in list(BPD_CACHE.items()):
        if not TL.has_edge(*k):
            del BPD_CACHE[k]
    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Cache to fold exterior loops.                                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
EFOLD_CACHE = {}

def fold_exterior_loop(seq, con, md, cache = None, spacer = 'NNN'):
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
        seq (str):          RNA sequence
        con (str):          RNA structure constraint
        md (RNA.md()):      ViennaRNA model details (temperature, noLP, etc.)
        spacer (str, optional): Specify the sequence of the spacer region.
                            Defaults to 'NNN'.

    Returns:
      (str, str):
    """
    if cache is None:
        global EFOLD_CACHE
    else:
        EFOLD_CACHE = cache


    pt = make_pair_table(con, base = 0)
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
    if ext_seq in EFOLD_CACHE:
        css = EFOLD_CACHE[ext_seq]
    else:
        fc_tmp = RNA.fold_compound(ext_seq, md)
        fc_tmp.constraints_add(ext_con, 
                RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP)
        css, cfe = fc_tmp.mfe()
        EFOLD_CACHE[ext_seq] = css
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

    return con, ext_seq

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
class TrafoLandscape(PrimePathLandscape):
    def __init__(self, *kargs, **kwargs):
        super(TrafoLandscape, self).__init__(*kargs, **kwargs)

        # Calculate full backtracking matrix for all subsequent MFE runs.
        if self.fc:
            _ = self.fc.mfe()

        # Private instance variables:
        self._transcript_length = 0
        self.total_time = 0

        # Default parameters:
        self.k0 = 2e5 # set directly
        self.minh = 0 # set using t_fast

    @property
    def transcript(self):
        return self.sequence[0:self._transcript_length]

    @property
    def t_fast(self):
        """Time-scale separation parameter for coarse-graining.

        Expected folding times faster than t-fast are considered instantaneous,
        and used to assign a conformation to a macrostate. (A folding time equal
        t-fast is not instantaneous, it corresponds to the minimum height of an
        energy barrier minh.) This parameter retruns a time scale based on the
        internal minh parameter, which quantifies the minimal energy barrier
        height to separate two macrostates.

        t_fast = 1/(self.k0 * exp(-self.minh/self.RT))
        """
        return 1/(self.k0 * math.exp(-self.minh/self.RT))

    @t_fast.setter
    def t_fast(self, value):
        self.minh = max(0, -self.RT * math.log(1 / value / self.k0))

    def print_structures(self):
        for ss in self.sorted_nodes():
            print(ss, self.nodes[ss]['energy'], 
                    self.nodes[ss]['identity'], 
                    self.nodes[ss]['active'], 
                    '{:.4f}'.format(self.nodes[ss]['occupancy']))
                    #self.nodes[ss]['hiddennodes'],
                    #self.nodes[ss]['lminreps'])

    def print_transitions(self, nodes = None):
        if nodes is None:
            nodes = self.nodes
        tl = len(self.transcript)

        for ss1, ss2 in combinations(nodes, 2):
            for [start, stop, fwb, rvb] in show_flooded_prime_path(self.sequence, ss1, ss2):
                print(start[:tl], stop[:tl], fwb, rvb)
            print()

    @property
    def active_subgraph(self):
        """ Return a new object with only active nodes (and edges) """
        other = super(TrafoLandscape, self).active_subgraph
        # Private instance variables:
        other._transcript_length = self._transcript_length
        other.total_time = self.total_time
        return other
 
    def expand(self, extend = 1, mfree = 6):
        """ Find new secondary structures and add them to :obj:`TrafoLandscape()`

        The function adds to types of new structures: 
            1) The mfe structure for the current sequence length.
            2) The helix-fraying of all currently active structures.

        Args:
            extend (int, optional): number of nucleotide extensions before graph
                expansion (updates the global variable transcript length). Defaults to 1.
            mfree (int, optional): minimum number of freed bases during a
                helix-opening step. Defaults to 6.

        Returns:
            int: Number of new nodes
        """
        rlog.info(f'Expand: {len(self.active_nodes)} extend = {extend} mfree = {mfree}')
        fseq = self.sequence
        self._transcript_length += extend
        if self._transcript_length > len(fseq):
            self._transcript_length = len(fseq)
        seq = self.transcript

        cnid = self.nodeID
        fc = self.fc

        # Calculate MFE of current transcript
        mfess, mfe = fc.backtrack(len(seq))
        future = '.' * (len(fseq) - len(seq))
        mfess = mfess + future

        # If there is no node because we are in the beginning, add the node.
        if len(self) == 0: 
            self.addnode(mfess, structure = mfess, active = True)
            self.nodes[mfess]['occupancy'] = 1
            nn = set([mfess])
        else: 
            # Save the set of active parent nodes for later.
            parents = self.active_nodes
            #nn = self.connect_nodes_n2(parents) 
            #rlog.info(f'Connected parents, new structs: {len(nn)}')
            #[rlog.debug(f' {new}') for new in nn]
            nn = self.connect_nodes_nm(nodesA = parents, 
                                       nodesB = [mfess], direct = False)
            rlog.info(f'Connected MFE: {mfess[:len(seq)]}, new structs: {len(nn)}')
            [rlog.info(f' {new}') for new in nn]

            # They can also be previously inactive nodes or other parents ...
            fraying_nodes = self.find_fraying_neighbors(parents, mfree = 6)
            for p, fn in fraying_nodes.items():
                # Connects fraying neighbors to other active nodes
                nn |= self.connect_nodes_nm(nodesA = self.active_nodes, nodesB = list(fn))
                # Connect fraying neighbors with each other
                nn |= self.connect_nodes_n2(nodes = list(fn), direct = True)
            rlog.info(f'Connected fraying neighbors, new structs: {len(nn)}')
            [rlog.info(f' {new}') for new in nn]
        return nn

    def find_fraying_neighbors(self, parents, mfree = 6):
        """ Test me!
        Returns:
            dict: new_nodes[parent]=set(neigbors)
        """
        seq = self.transcript
        future = '.' * (len(self.sequence) - len(seq))
        md = self.md
        # Keep track of all nodes and how they are related to parents.
        new_nodes = dict() 
        for ni in parents:
            # short secondary structure (without its future)
            sss = ni[0:len(seq)]
            # compute a set of all helix fraying open steps
            opened = open_fraying_helices(seq, sss, mfree)
            # Do a constrained exterior loop folding for all fraying structures
            connected = set([ni]) # add the parent
            for onbr in opened:
                nbr, _ = fold_exterior_loop(seq, onbr, md)
                nbr += future
                connected.add(nbr)
            new_nodes[ni] = connected
        return new_nodes

    def get_simulation_files_tkn(self, basename, snodes = None):
        """ Print a rate matrix and the initial occupancy vector.

        Overwrites: Ribolands.get_simulation_files_tkn()

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

        """
        seq = self.transcript
        if snodes is None:
            snodes = self.sorted_nodes(attribute = 'energy') 

        num = len(snodes)+1

        bofile = basename + '_lands.bar'
        brfile = basename + '_rates.txt'
        bbfile = basename + '_rates.bin'
        p0 = []

        with open(bofile, 'w') as bar, open(brfile, 'w') as rts, open(bbfile, 'wb') as brts:
            bar.write("  ID {}  Energy  {}\n".format(seq, 
                ' '.join(map("{:7d}".format, range(1, num)))))
            brts.write(pack("i", len(snodes)))
            for ni, node in enumerate(snodes, 1):
                ne = self.nodes[node]['energy']
                no = self.nodes[node]['occupancy']

                # Calculate barrier heights to all other basins.
                barstr = ''
                for other in snodes:
                    oe = self.nodes[other]['energy']
                    sE = self.get_saddle(node, other)
                    if sE is not None:
                        barstr += ' {:7.2f}'.format(sE - ne)
                    else:
                        barstr += ' {:7.2f}'.format(float('nan'))

                # Print structures and neighbors to bfile:
                bar.write("{:4d} {} {:7.2f} {}\n".format(ni, node[:len(seq)], ne, barstr))

                # Add ni occupancy to p0
                if no > 0:
                    p0.append("{}={}".format(ni, no))

                # Print rate matrix to rfile and brfile
                trates = []
                rates = []
                for other in snodes:
                    if self.has_edge(node, other):
                        rates.append(self[node][other]['weight'])
                        trates.append(self[other][node]['weight'])
                    else:
                        rates.append(0)
                        trates.append(0)

                line = "".join(map("{:10.4g}".format, rates))
                rts.write("{}\n".format(line))
                for r in trates:
                    brts.write(pack("d", r))

        return bbfile, brfile, bofile, p0

    def update_occupancies_tkn(self, tfile, snodes):
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
                    ss = snodes[e]
                    self.nodes[ss]['occupancy'] = float(occu)/tot_occ
            else :
                time = float(lastlines[0].split()[0])
                iterations = int(lastlines[-1].split()[-1])
                tot_occ = sum(map(float, lastlines[0].split()[1:]))
                for e, occu in enumerate(lastlines[0].split()[1:]):
                    ss = snodes[e]
                    assert self.nodes[ss]['active']
                    sso = self.nodes[ss]['occupancy']
                    self.nodes[ss]['occupancy'] = float(occu)/tot_occ
        return time, iterations

    def prune(self, pmin, keepers = None):
        """ 
        We prune by base-pair distance to active nodes... 

        Identify the nodes we have to keep, either because they are
        occupied local minima, or because they are basins connecting
        occupied local minima.

        """
        rlog.info(f'Pruning: pmin = {pmin} keepers = {len(keepers) if keepers else 0}')
        if keepers is None:
            keepers = set()
            for n in self.active_nodes:
                no = self.nodes[n]['occupancy']
                if no >= pmin:
                    keepers.add(n)
                    rlog.debug(f'Keeper: {n}')

        if len(keepers) == 1:
            rlog.info('Only one active node with occupancy!')
            for node in list(self.nodes):
                self.nodes[node]['pruned'] = self.nodes.get('pruned', 0)
                if node in keepers:
                    assert self.nodes[node]['lminreps'] is None
                    self.nodes[node]['active'] = True
                    self.nodes[node]['pruned'] = 0
                else:
                    self.nodes[node]['active'] = False
                    self.nodes[node]['pruned'] += 1
        else:
            # Inactivate old nodes and increment "pruned"
            for x, y in combinations(keepers, 2):
                if self.has_edge(x, y):
                    continue
                for [start, stop, fwb, fwg, rvb] in self.get_prime_path_minima(x, y):
                    if stop not in self.nodes:
                        rlog.debug(f'Minimal prime path basins are not directly connected!\n' + \
                                f'{x}\n{y}')
                        break
                    d1 = self.nodes[start]
                    d2 = self.nodes[stop]
                    basin1 = set([start]) if d1['lminreps'] is None else d1['lminreps']
                    basin2 = set([stop]) if d2['lminreps'] is None else d2['lminreps']
                    if not any(s1 == s2 or self.get_barrier(s1, s2) is not None \
                            for (s1, s2) in product(basin1, basin2)):
                        rlog.debug(f'Minimal prime path basins are not directly connected!\n' + \
                                f'{x}\n{y}')
                    keepers |= basin1 | basin2
            for node in list(self.nodes):
                self.nodes[node]['pruned'] = self.nodes.get('pruned', 0)
                if node in keepers:
                    assert self.nodes[node]['lminreps'] is None
                    self.nodes[node]['active'] = True
                    self.nodes[node]['pruned'] = 0
                else:
                    self.nodes[node]['active'] = False
                    self.nodes[node]['pruned'] += 1

        parents = self.active_nodes
        for node in self.sorted_nodes(attribute = 'energy', rev = True):
            self.nodes[node]['lminreps'] = None
            self.nodes[node]['hiddennodes'] = None
            if not self.nodes[node]['active'] and self.nodes[node]['occupancy'] > 0:
                assert self.nodes[node]['pruned'] == 1
                # We prune by base-pair distance to active nodes... 
                move = sorted(parents, key = lambda x: get_bpd_cache(node, x))
                self.nodes[move[0]]['occupancy'] += self.nodes[node]['occupancy']
                self.nodes[node]['occupancy'] = 0
            if self.nodes[node]['pruned'] > 0:
                rlog.debug(f'Removing node: {node}')
                self.remove_node(node)
        return

    def sorted_trajectories_iter(self, snodes, tfile, softmap = None):
        """ Yields the time course information using a treekin output file.

        Args:
          snodes (list): a list of nodes sorted by their energy
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
                            ss = snodes[e][0]
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
                        ss = snodes[e][0]
                        sss = ss[0:self._transcript_length]
                        yield self.nodes[ss]['identity'], ttime + time, occu, \
                            sss, self.nodes[ss]['energy']
                    prevcourse = course
        return

def open_fraying_helices(seq, ss, free = 6):
    """ Generate structures with opened fraying helices. 
    
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
    pt = make_pair_table(ss, base = 0)

    # mutable secondary structure
    nbr = list(ss)

    def rec_fill_nbrs(ss, mb, pt, myrange, free):
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
                    rec_fill_nbrs(''.join(nb), mb, pt, (p, q), free - o)
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

    rec_fill_nbrs(ss, nbr, pt, (0, len(ss)), free)
    nbrs.add(''.join(nbr))
    return nbrs


