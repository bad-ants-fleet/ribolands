#
# ribolands.trafo
# 
# Home of the TrafoLandscape object.
#

import logging
rlog = logging.getLogger(__name__)

import os
import sys
import math
from struct import pack
from datetime import datetime
from itertools import combinations, product, islice

import RNA
from ribolands import RiboLandscape
from ribolands.utils import make_pair_table
from ribolands.pathfinder import (BPD_CACHE, get_bpd_cache, 
                                  get_guide_graph,
                                  neighborhood_flooding,
                                  top_down_coarse_graining)

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

def find_fraying_neighbors(seq, md, parents, mfree = 6):
    """ 
    Returns:
        dict: new_nodes[parent]=set(neigbors)
    """
    future = '.' * (len(parents[0]) - len(seq))
    # Keep track of all nodes and how they are related to parents.
    new_nodes = dict() 
    for ni in parents:
        # short secondary structure (without its future)
        sss = ni[0:len(seq)]
        assert len(sss + future) == len(parents[0])
        # compute a set of all helix fraying open steps
        opened = open_fraying_helices(seq, sss, mfree)
        # Do a constrained exterior loop folding for all fraying structures
        connected = set() 
        for onbr in opened:
            nbr, _ = fold_exterior_loop(seq, onbr, md)
            nbr += future
            connected.add(nbr)
        connected.discard(ni) # remove the parent
        new_nodes[ni] = connected
    return new_nodes


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
class TrafoLandscape(RiboLandscape):
    def __init__(self, *kargs, **kwargs):
        super(TrafoLandscape, self).__init__(*kargs, **kwargs)

        # Calculate full backtracking matrix for all subsequent MFE runs.
        _ = self.fc.mfe()
        #self.fp = init_findpath_max(self.sequence)

        # Private instance variables:
        self._transcript_length = 0
        self._cg_edges = dict()

        # TODO
        self.total_time = 0

        # Default parameters:
        self.prefix = '' # for autogenerated node IDs
        self.k0 = 2e5 # set directly
        self.minh = 0 # [dcal/mol] set using t_fast
        self.maxh = 0 # [dcal/mol] set using t_slow
        self.fpath = 0

    @property
    def transcript(self):
        return self.sequence[0:self._transcript_length]

    def addnode(self, *kargs, **kwargs):
        """ Add a node with specific tags:

        - Active (bool) is to filter a subset of interesting nodes.
        - Lminreps (set) is to return a set of nodes that are the local minimum
            reperesentatives of this node.
        - Hiddennodes (set) is to return a set of nodes that are associated with
            this local minimum.
        After coarse graining, a node should have either lminreps or
        hiddennodes set, but not both. 
        """
        if 'active' not in kwargs:
            kwargs['active'] = True
        if 'pruned' not in kwargs:
            kwargs['pruned'] = 0
        if 'lminreps' not in kwargs:
            kwargs['lminreps'] = None
        if 'hiddennodes' not in kwargs:
            kwargs['hiddennodes'] = None
        if 'occtransfer' not in kwargs:
            kwargs['occtransfer'] = None
        super(TrafoLandscape, self).addnode(*kargs, **kwargs)

    @property
    def cg_edges(self):
        return self._cg_edges

    def has_edge(self, s1, s2, cg = False):
        # Function overload to allow a set of coarse grained edges.
        return ((s1, s2) in self.cg_edges) if cg else ((s1, s2) in self.edges)

    def get_saddle(self, s1, s2):
        """ Returns the saddle energy of a transition edge.  """
        return self.edges[(s1, s2)]['saddle_energy'] if self.has_edge(s1, s2, cg = False) else None

    def get_cg_saddle(self, s1, s2):
        """ Returns the saddle energy of a coarse grained transition edge.  """
        return self.cg_edges[(s1, s2)]['saddle_energy'] if self.has_edge(s1, s2, cg = True) else None

    @property
    def local_mins(self):
        for n in self.nodes:
            if not self.nodes[n]['lminreps']:
                yield n

    @property
    def active_local_mins(self):
        for n in self.nodes:
            if self.nodes[n]['active'] and not self.nodes[n]['lminreps']:
                yield n
    @property
    def hidden_nodes(self):
        for n in self.nodes:
            if self.nodes[n]['lminreps']:
                yield n

    @property
    def active_nodes(self):
        for n in self.nodes:
            if self.nodes[n]['active']:
                yield n

    @property
    def inactive_nodes(self):
        for n in self.nodes:
            if not self.nodes[n]['active']:
                assert self.nodes[n]['active'] is False
                yield n

    def expand(self, extend = 1, mfree = 6):
        """ Find new secondary structures and determine their neighborhood in the landscape.

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
        fseq = self.sequence
        self._transcript_length += extend
        if self._transcript_length > len(fseq):
            self._transcript_length = len(fseq)
        seq = self.transcript
        fc = self.fc

        # Calculate MFE of current transcript.
        mfess, _ = fc.backtrack(len(seq))
        future = '.' * (len(fseq) - len(seq))
        mfess = mfess + future

        i_time = datetime.now()

        # If there is no node because we are in the beginning, add the node.
        if len(self.nodes) == 0: 
            self.addnode(mfess, structure = mfess, occupancy = 1)
            nn = set([mfess])
            on = set()
        else: 
            md = self.md
            parents = list(self.active_local_mins)
            fraying_nodes = find_fraying_neighbors(seq, md, parents, mfree = 6)

            # 1) Add all new structures to the set of nodes.
            nn, on = set(), set()
            if mfess not in self.nodes:
                self.addnode(mfess, structure = mfess)
                nn.add(mfess)
            elif not self.nodes[mfess]['active']:
                on.add(mfess)
            self.nodes[mfess]['active'] = True
            for fns in fraying_nodes.values():
                for fn in fns:
                    if fn not in self.nodes:
                        self.addnode(fn, structure = fn)
                        nn.add(fn)
                    elif not self.nodes[fn]['active']:
                        on.add(fn)
                    self.nodes[fn]['active'] = True

            f_time = datetime.now()

            ndata = {n[0:len(seq)]: d for n, d in self.nodes.items() if d['active']} 
            gnodes, gedges = get_guide_graph(seq, md, ndata.keys())
            assert all(ss != '' for (ss, en) in gnodes)

            for (ss, en) in gnodes:
                if ss+future not in self.nodes:
                    self.addnode(ss+future, structure = ss+future)
                    nn.add(ss+future)
                elif not self.nodes[ss+future]['active']:
                    on.add(ss+future)
                self.nodes[ss+future]['active'] = True
            ndata = {n: d for n, d in self.nodes.items() if d['active']} 

            lgedges = set()
            for (x, y) in gedges:
                lgedges.add((x+future, y+future))
            gedges = lgedges

            g_time = datetime.now()

            # 2) Include edge-data from previous network if nodes are active.
            edata = {k: v for k, v in self.edges.items() if (self.nodes[k[0]]['active'] and self.nodes[k[1]]['active']) and 
                    v['saddle_energy'] is not None}

            ndata, edata = neighborhood_flooding((fseq, md), ndata, gedges, tedges = edata, minh = self.minh)
            #ndata, edata = neighborhood_flooding(self.fp, ndata, gedges, tedges = edata, minh = self.minh)

            # 3) Extract new node data.
            for node in ndata:
                if node not in self.nodes:
                    self.addnode(node, structure = node, energy = ndata[node]['energy'])
                    nn.add(node)
                elif not self.nodes[node]['active']:
                    on.add(node)
                self.nodes[node]['active'] = True
                assert self.nodes[node]['energy'] == ndata[node]['energy']

            # 4) Update to new edges.
            for (x, y) in edata:
                se = edata[(x, y)]['saddle_energy']
                self.addedge(x, y, saddle_energy = se)

            l_time = datetime.now()
            frayytime = (f_time - i_time).total_seconds() 
            guidetime = (g_time - f_time).total_seconds()
            floodtime = (l_time - g_time).total_seconds()
            tottime = (l_time - i_time).total_seconds()
            #print(len(seq), tottime, frayytime, guidetime, floodtime)
            #sys.stdout.flush()
        return nn, on

    def get_coarse_network(self, minh = None):
        """ Produce a smaller graph of local minima and best connections.

        active vs inactive: is useful because we might have to assign IDs for simulation.
        Structures disappear and reappear, so we keep a history of inactive nodes before
        finally deleting them.
        """
        if minh is None:
            minh = self.minh

        ndata = dict()
        for n in self.nodes: 
            self.nodes[n]['lminreps'] = set()
            self.nodes[n]['hiddennodes'] = set()
            # Because this node remained inactive during graph expansion, we
            # can now safely transfer its occupancy.
            if not self.nodes[n]['active'] and self.nodes[n]['occupancy'] != 0:
                for tn in self.nodes[n]['occtransfer']:
                    assert self.nodes[tn]['active']
                    self.nodes[tn]['occupancy'] += self.nodes[n]['occupancy']/len(self.nodes[n]['occtransfer'])
                self.nodes[n]['occupancy'] = 0
        ndata = {n: d for n, d in self.nodes.items() if d['active']} # only active.
        edata = {k: v for k, v in self.edges.items() if self.nodes[k[0]]['active'] and self.nodes[k[1]]['active']}
        cg_ndata, cg_edata, cg_mapping = top_down_coarse_graining(ndata, edata, minh)
        assert all((n in ndata) for n in cg_ndata)

        # Translate coarse grain results to TL.
        self._cg_edges = dict()
        for (x, y) in cg_edata:
            se = cg_edata[(x, y)]['saddle_energy']
            ex = cg_ndata[x]['energy']
            bar = (se-ex) / 100
            self._cg_edges[(x, y)] = {'saddle_energy': se,
                                      'weight': self.k0 * math.e**(-bar/self.RT)}

        for lmin, hidden in cg_mapping.items():
            assert self.nodes[lmin]['active']
            for hn in hidden:
                assert self.nodes[hn]['active']
                self.nodes[hn]['lminreps'].add(lmin)
            self.nodes[lmin]['hiddennodes'] = hidden

        # Move occupancy to lmins.
        for hn in self.hidden_nodes:
            if not self.nodes[n]['active']:
                assert self.nodes[n]['occupancy'] == 0
            if self.nodes[hn]['occupancy']:
                for lrep in self.nodes[hn]['lminreps']:
                    self.nodes[lrep]['occupancy'] += self.nodes[hn]['occupancy']/len(self.nodes[hn]['lminreps'])
            self.nodes[hn]['occupancy'] = 0
            # TODO: try this ...
            self.nodes[hn]['active'] = False
        return len(cg_ndata), len(cg_edata)

    def get_simulation_files_tkn(self, basename):
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
        snodes = sorted(self.active_local_mins, key = lambda x: self.nodes[x]['energy'])
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
                ne = self.nodes[node]['energy']
                no = self.nodes[node]['occupancy']

                # Calculate barrier heights to all other basins.
                barstr = ''
                for other in snodes:
                    oe = self.nodes[other]['energy']
                    sE = self.get_cg_saddle(node, other)
                    if sE is not None:
                        barstr += ' {:7.2f}'.format((sE - ne)/100)
                    else:
                        barstr += ' {:7.2f}'.format(float('nan'))

                # Print structures and neighbors to bfile:
                bar.write("{:4d} {} {:7.2f} {}\n".format(ni, node[:len(seq)], ne/100, barstr))

                # Add ni occupancy to p0
                if no > 0:
                    p0.append("{}={}".format(ni, no))

                # Print rate matrix to rfile and brfile
                trates = []
                rates = []
                for other in snodes:
                    if self.has_edge(node, other, cg = True):
                        rates.append(self.cg_edges[(node, other)]['weight'])
                        trates.append(self.cg_edges[(other, node)]['weight'])
                    else:
                        rates.append(0)
                        trates.append(0)
                line = "".join(map("{:10.4g}".format, rates))
                rts.write("{}\n".format(line))
                for r in trates:
                    brts.write(pack("d", r))
        return snodes, bbfile, brfile, bofile, p0

    def update_occupancies_tkn(self, tfile, snodes):
        """ Update the occupancy in the graph and return simulation time.
        """
        with open(tfile, 'rb') as f:
            f.seek(-2, os.SEEK_END)
            skipchar = False
            while True:
                if not skipchar and f.read(1) == b'\n':
                    fchar = f.read(1)
                    if fchar == b'#':
                        skipchar = True
                    else:
                        break
                else:
                    skipchar = False
                f.seek(-2, os.SEEK_CUR)
            last_line = fchar.decode() + f.readline().decode()

        ll = list(map(float, last_line.split()))
        time, occus = ll[0], ll[1:]
        tot_occ = sum(occus)
        for e, occu in enumerate(occus):
            ss = snodes[e]
            self.nodes[ss]['occupancy'] = occu/tot_occ
        return time

    def prune(self, pmin, delth = 3):
        """ TODO make sure return values make sense... distinguish lmins and hn.
        """
        new_inactive_lms = []
        tot_pruned = 0
        for lm in sorted(self.active_local_mins, key = lambda x: self.nodes[x]['occupancy']):
            self.nodes[lm]['pruned'] = 0
            if tot_pruned + self.nodes[lm]['occupancy'] > pmin:
                break
            tot_pruned += self.nodes[lm]['occupancy']
            for hn in self.nodes[lm]['hiddennodes']:
                assert self.nodes[hn]['occupancy'] == 0
                self.nodes[hn]['active'] = False
            self.nodes[lm]['active'] = False
            new_inactive_lms.append(lm)

        def get_active_nbrs(lm, forbidden = None):
            # TODO: lazy and probably inefficient implementation.
            if forbidden is None:
                forbidden = set()
            forbidden.add(lm)
            found = False
            remaining = []
            for (x, y) in self.cg_edges:
                assert x != y
                if x == lm and y not in forbidden:
                    if self.nodes[y]['active']:
                        found = True
                        yield y
                    else:
                        remaining.append(y)
            if not found:
                for r in remaining:
                    assert self.nodes[r]['active'] is False
                    for a in get_active_nbrs(r, forbidden):
                        yield a
            return

        for lm in new_inactive_lms:
            assert not self.nodes[lm]['active']
            # TODO: should it really be all of them?
            self.nodes[lm]['occtransfer'] = set(get_active_nbrs(lm))
            #print('occu', lm, '->',  self.nodes[lm]['occtransfer'])
            assert len(self.nodes[lm]['occtransfer']) > 0

        pn = set()
        dn = set()
        for node in list(self.nodes):
            if self.nodes[node]['active']:
                self.nodes[node]['pruned'] = 0
            else:
                self.nodes[node]['pruned'] += 1
                if self.nodes[node]['pruned'] == 1:
                    # this includes hidden nodes that were not part of the simulation
                    pn.add(node)
            rlog.debug(f'After pruning: {node} {self.nodes[node]}')
            if self.nodes[node]['pruned'] > delth:
                #TODO: probably quite inefficient...
                for (x, y) in set(self.edges):
                    if node in (x, y):
                        del self._edges[(x,y)]
                del self._nodes[node]
                dn.add(node)
        return pn, dn

