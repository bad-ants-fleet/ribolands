#!/usr/bin/env python
#
# ribolands.pathfinder
# 
# Utilities using the findpath algorithm.
#

import logging
rlog = logging.getLogger(__name__)

import RNA
import math
import numpy as np
from itertools import chain, combinations, product

from .utils import make_pair_table

MAXPATH = True
try: # Check if findpath python module is installed
    from .maxpath import init_findpath_max, findpath_max
except ImportError:
    MAXPATH = False

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Custom error definitions                                                     #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
class PathfinderError(Exception):
    pass

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Cache to look-up base-pair distances                                         #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
BPD_CACHE = {}

def clear_bpd_cache():
    global BPD_CACHE
    BPD_CACHE = {}

def get_bpd_cache(s1, s2):
    global BPD_CACHE
    if (s1, s2) in BPD_CACHE:
        return BPD_CACHE[(s1, s2)]
    elif (s2, s1) in BPD_CACHE:
        return BPD_CACHE[(s2, s1)]
    else:
        dist = RNA.bp_distance(s1, s2)
        BPD_CACHE[(s1, s2)] = dist
        return dist

BPD_I_CACHE = {}
def get_bpd_i_cache(p, q):
    global BPD_I_CACHE
    def intersect(p, q):
        ptp = RNA.ptable(p)
        ptq = RNA.ptable(q)
        bpp = set((i, j) for i, j in enumerate(ptp[1:], 1) if j > i)
        bpq = set((i, j) for i, j in enumerate(ptq[1:], 1) if j > i)
        bpI = bpp & bpq # intersect has always less basepairs
        dpI = len(bpp) - len(bpI)
        dIq = len(bpq) - len(bpI)
        return dpI, dIq
    if (p, q) not in BPD_I_CACHE:
        assert (q, p) not in BPD_I_CACHE
        dpI, dIq = intersect(p, q)
        BPD_I_CACHE[(p, q)] = dpI
        BPD_I_CACHE[(q, p)] = dIq
    return BPD_I_CACHE[(p, q)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Findpath stuff                                                               #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def common_exterior_bases(pt1, pt2):
    hide = -1 
    for e in range(0, len(pt1)):
        if e > hide and pt1[e] == -1 and pt2[e] == -1: # a common unpaired base
            yield e
        elif pt1[e] > e or pt2[e] > e: # a common base-pair
            hide = max([hide, pt1[e], pt2[e]])

def common_basepairs(pt1, pt2):
    return ((i, j1) for i, (j1, j2) in enumerate(zip(pt1, pt2)) if j1 == j2 and i < j1)

def split_struct(struct, i, j, spacer = '...'):
    if j is None:
        return struct[:i+1], struct[i:]
    else:
        return struct[:i+1] + spacer + struct[j:], struct[i:j+1]

def merge_struct(outside, inside, i, j, slen = 4):
    if j is None:
        return outside[:i] + inside
    else:
        return outside[:i] + inside[:-1] + outside[j-len(inside)+slen+1:]

def findpath_merge(outside, inside, i, j):
    """ Merge two composable findpath runs.

    This is a simplified variant which just appends the two runs ...
    """
    starten = outside[0][1] + inside[0][1]
    bh1 = bh2 = 0

    path1 = []
    stepI, enI = inside[0]
    for stepO, enO in outside:
        merged = merge_struct(stepO, stepI, i, j)
        energy = enO + enI - starten
        path1.append((merged, energy))
        bh1 = max(bh1, energy)
    stepO, enO = outside[-1]
    for stepI, enI in inside[1:]:
        merged = merge_struct(stepO, stepI, i, j)
        energy = enO + enI - starten
        path1.append((merged, energy))
        bh1 = max(bh1, energy)

    path2 = []
    stepO, enO = outside[0]
    for stepI, enI in inside:
        merged = merge_struct(stepO, stepI, i, j)
        energy = enO + enI - starten
        path2.append((merged, energy))
        bh2 = max(bh2, energy)
    stepI, enI = inside[-1]
    for stepO, enO in outside[1:]:
        merged = merge_struct(stepO, stepI, i, j)
        energy = enO + enI - starten
        bh2 = max(bh2, energy)
        path2.append((merged, energy))
    return (path1, bh1) if bh1 < bh2 else (path2, bh2)

def findpath_split(seq, ss1, ss2, md, th = 5, w = None):
    """ Calculate findpath barriers for smaller components.

    Args:
        seq: RNA sequence.
        ss1: Structure 1.
        ss2: Structure 2.
        md: ViennaRNA model details.
        th: Threshold of how many basepairs must change for an independent findpath run. Defaults to 5.
        w: Findpath width. Defaults to None.

    Returns:
        path, barrier: The folding path and the barrier height. 
            WARNING: If path splitting actually took place, then energy values
            given in the path data are only relative to the starting structure.
    """
    pt1 = make_pair_table(ss1, base = 0, chars = list('.x'))
    pt2 = make_pair_table(ss2, base = 0, chars = list('.x'))
    mindiff = None
    recurse = None
    for ij in chain(common_exterior_bases(pt1, pt2),
                    common_basepairs(pt1, pt2)):
        (i, j) = ij if isinstance(ij, tuple) else (ij, None)
        st1O, st1I = split_struct(ss1, i, j, spacer = '...')
        st2O, st2I = split_struct(ss2, i, j, spacer = '...')
        do = RNA.bp_distance(st1O, st2O)
        if do < th: continue
        di = RNA.bp_distance(st1I, st2I)
        if di < th: continue
        diff = abs(di-do)
        if mindiff is None or diff < mindiff:
            mindiff = diff
            seqO, seqI = split_struct(seq, i, j, spacer = 'NNN')
            recurse = ((i, j),
                       (seqO, st1O, st2O), 
                       (seqI, st1I, st2I))
        elif mindiff is not None and diff > mindiff:
            # No need to check the rest if we are getting worse.
            break

    if mindiff is not None:
        pathO, _ = findpath_split(*recurse[1], md, th, w)
        pathI, _ = findpath_split(*recurse[2], md, th, w)
        return findpath_merge(pathO, pathI, *recurse[0])
    else:
        fpw = 4 * RNA.bp_distance(ss1, ss2) if w is None else w
        return call_findpath(seq, ss1, ss2, md, w = fpw)

def call_findpath(seq, ss1, ss2, md, w, maxbar = float('inf')):
    """ Call ViennaRNA findpath.

    TODO: This is where we could minimize the motif and look it up in a cache.
    """
    fc = RNA.fold_compound(seq, md)
    if maxbar == float('inf'):
        path = fc.path_findpath(ss1, ss2, width = w)
    else:
        e1 = round(fc.eval_structure(ss1), 2)
        dcal_bound = int(round((maxbar + e1) * 100))
        path = fc.path_findpath(ss1, ss2, maxE = dcal_bound, width = w)
    del fc

    if len(path):
        mypath = []
        barrier = None
        for step in path:
            struct = step.s
            energy = int(round(step.en*100))
            mypath.append((struct, energy))
            if barrier is None or barrier < energy:
                barrier = energy
        barrier -= int(round(path[0].en*100))
        del step, path # potentially a good idea.
        return mypath, barrier
    return None, None

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Path flooding                                                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def path_flooding(path, minh, maxlm = None):
    """ Use flooding algorithm to determine local minima on a folding path.
    
    Identifies the lowest energy of the local minimum, one representative
    stucture and "other structures" associated with that minimum. Beware that
    the "other structures" can contain saddle components when the function
    is applied to paths with degenerate saddle components.

    Args:
        path (list[(str, int)]): A list of tuples (structure, energy). 
            The energy must be an integer with unit dcal/mol.
        minh (int): Minimal height of an energy barrier separatig two
            basins [dcal/mol].
        maxlm (optional, int): The highest possible energy of a local minimum.

    Returns:
        [dict, dict]: properties of the local minima.
    """
    assert isinstance(minh, int) # dcal/mol
    assert isinstance(path[0][1], int) # dcal/mol
    indexpath = [[ss, en, step] for step, [ss, en] in enumerate(path)]
    sconf = sorted(indexpath, key = lambda x: x[1]) # by energy

    lmins = dict() # lmins[ss] = [en, set(steps)]
    ssmap = dict() # ssmap[step] = step
    for e, (ss, en, step) in enumerate(sconf):
        if step == 0 or step == len(path) - 1: # edge cases
            nbr = step + 1 if step == 0 else step - 1
            if nbr in ssmap:
                nbr = ssmap[nbr]
                lms, lme = path[nbr]
                assert en >= lme
                assert lme == lmins[lms][0]
                lmins[lms][1].add(step)
                ssmap[step] = nbr
            else:
                lmins[ss] = [en, set([step])]
                ssmap[step] = step
        else: 
            assert 0 < step < len(path)-1
            left, right = step-1, step+1
            if left in ssmap and right in ssmap: # let's connect!
                left, right = ssmap[left], ssmap[right]
                assert left != right # yay, saddle!
                lmsl, lmel = path[left]
                lmsr, lmer = path[right]
                [lower, higher] = [left, right] if lmel < lmer else [right, left]
                [lowlm, highlm] = [lmsl, lmsr] if lmel < lmer else [lmsr, lmsl]
                [lowen, highen] = [lmel, lmer] if lmel < lmer else [lmer, lmel]
                if en - highen < minh or (maxlm and highen > maxlm): # merge higher minimum into lower one.
                    lmins[lowlm][1] |= lmins[highlm][1] | set([step]) 
                    for hlm in lmins[highlm][1]:
                        ssmap[hlm] = lower
                    ssmap[step] = lower
                    del lmins[highlm]
                else: # keep them separated ...
                    ssmap[step] = [higher, lower]
            elif left in ssmap:
                lms, lme = path[ssmap[left]]
                assert en >= lme
                lmins[lms][1].add(step)
                ssmap[step] = ssmap[left]
            elif right in ssmap:
                lms, lme = path[ssmap[right]]
                assert en >= lme
                lmins[lms][1].add(step)
                ssmap[step] = ssmap[right]
            else:
                lmins[ss] = [en, set([step])]
                ssmap[step] = step
    return ssmap

def edge_flooding(fp, s1, s2, e1, e2, minh = None):
    """ Connect two arbitrary secondary structures.
        Args:
            minh (int): Minimal height of an energy barrier separatig 
                        two basins in dcal/mol.
    """
    if isinstance(fp, tuple):
        (seq, md) = fp
        path, barrier = findpath_split(seq, s1, s2, md)
        if path[0][1] == 0:
            for i in range(len(path)):
                path[i] = (path[i][0], path[i][1] + e1)
        sen = e1 + barrier
    else:
        sen, path = findpath_max(fp, s1, s2)

    if path[0][1] != e1:
        assert path[0][1] == e1

    if minh is not None:
        # NOTE: We only allow to add minima with energy between start/stop structure.
        maxlme = max(path[0][1], path[-1][1])
        ssmap = path_flooding(path, minh, maxlm = maxlme)
        for e, si in enumerate(sorted(ssmap)):
            lm = ssmap[si]
            if isinstance(lm, list):
                [si1, si2] = sorted(lm)
                (ssB, enB) = path[si]
                (ss1, en1) = path[si1]
                (ss2, en2) = path[si2]
                yield ss1, en1, ssB, enB, ss2, en2
            elif e == 0 and si != lm:
                #rlog.warning(f'Starting structure {s1} is not a mimimum!')
                (ss2, en2) = path[lm]
                yield s1, e1, s1, e1, ss2, en2
            elif e == len(ssmap) - 1 and si != lm:
                #rlog.warning(f'Final structure {s2} is not a mimimum!')
                (ss1, en1) = path[lm]
                yield ss1, en1, s2, e2, s2, e2
    else:
        yield s1, e1, None, sen, s2, e2

def neighborhood_flooding(fp, ndata, gedges, tedges = None, minh = None):
    """ Calculate flooded paths for all edges.

    Note: Modifies arguments ndata and tedges.

    Starting with ndata and tedges (which may be None), 
    we add all the guide-edges to the graph.

    Args:
        ndata:
        gedges:
        tedges: 
    """
    if tedges is None:
        tedges = dict()

    guide_nbrs = {k: set() for k in ndata}
    for (p, q) in gedges: 
        assert p != q
        if (p, q) not in tedges:
            guide_nbrs[p].add(q)

    tstep_nbrs = {k: set() for k in ndata}
    for (p, q) in tedges: 
        assert p != q
        assert tedges[(p, q)]['saddle_energy'] is not None
        tstep_nbrs[p].add(q)

    while not all(len(v) == 0 for v in guide_nbrs.values()):
        seen = set()
        new_gedges = []
        for s2 in sorted(ndata, key = lambda x: (ndata[x]['energy'], x), 
                         reverse = True):
            for s1 in sorted(guide_nbrs[s2], key = lambda x: (ndata[x]['energy'], x), 
                             reverse = True):
                assert s1 != s2
                if (s1, s2) in seen:
                    continue
                assert ndata[s1]['energy'] <= ndata[s2]['energy']
                e2 = ndata[s2]['energy']
                e1 = ndata[s1]['energy']
                for (ss2, en2, ssB, enB, ss1, en1) in edge_flooding(fp, s2, s1, 
                                                                    e2, e1, minh):
                    assert ss1 != ss2
                    if ss2 == s2 and ss1 == s1: 
                        # adding the direct connection.
                        assert ndata[ss2]['energy'] == en2
                        assert ndata[ss1]['energy'] == en1
                        tedges[(ss2, ss1)] = {'saddle_energy': enB}
                        tedges[(ss1, ss2)] = {'saddle_energy': enB}
                        tstep_nbrs[ss2].add(ss1)
                        tstep_nbrs[ss1].add(ss2)
                        guide_nbrs[ss2].remove(ss1)
                        guide_nbrs[ss1].remove(ss2)
                    elif ss2 in tstep_nbrs.get(ss1, set()): 
                        # we might have found a better transition energy.
                        assert ss1 in tstep_nbrs[ss2]
                        # TODO: is this ncecessary?
                        #assert ss2 not in guide_nbrs[ss1]
                        #assert ss1 not in guide_nbrs[ss2]
                        assert ndata[ss2]['energy'] == en2
                        assert ndata[ss1]['energy'] == en1
                        if tedges[(ss2, ss1)]['saddle_energy'] is not None:
                            enB = min(enB, tedges[(ss2, ss1)]['saddle_energy'])
                        tedges[(ss2, ss1)] = {'saddle_energy': enB}
                        tedges[(ss1, ss2)] = {'saddle_energy': enB}
                    else: 
                        # postpone evaluation for next time.
                        assert ss2 not in ndata or ndata[ss2]['energy'] == en2
                        assert ss1 not in ndata or ndata[ss1]['energy'] == en1
                        ndata[ss2] = ndata.get(ss2, {'energy': en2})
                        ndata[ss1] = ndata.get(ss1, {'energy': en1})
                        guide_nbrs[ss2] = guide_nbrs.get(ss2, set())
                        guide_nbrs[ss1] = guide_nbrs.get(ss1, set())
                        tstep_nbrs[ss2] = tstep_nbrs.get(ss2, set())
                        tstep_nbrs[ss1] = tstep_nbrs.get(ss1, set())
                        new_gedges.append((ss2, ss1))
                if s2 not in tstep_nbrs[s1]: 
                    assert s1 != s2
                    assert s1 not in tstep_nbrs[s2] 
                    guide_nbrs[s2].remove(s1)
                    guide_nbrs[s1].remove(s2)
                    tedges[(s1, s2)] = {'saddle_energy': None}
                    tedges[(s2, s1)] = {'saddle_energy': None}
                seen.add((s2, s1))
        for (p, q) in new_gedges:
            assert p != q
            guide_nbrs[p].add(q)
            guide_nbrs[q].add(p)
    return ndata, tedges

def top_down_coarse_graining(ndata, edata, minh = 0):
    cg_ndata = {k: v for (k, v) in ndata.items()}
    cg_edata = {k: v for (k, v) in edata.items() if v['saddle_energy'] is not None}
    cg_basin = dict()

    successors = {k: set() for k in ndata}
    for (x, y) in cg_edata:
        assert (y, x) in cg_edata
        successors[x].add(y)

    for node in sorted(ndata, key = lambda x: (ndata[x]['energy'], x), reverse = True):
        en = cg_ndata[node]['energy']
        nbrs = successors[node]
        for nb in nbrs: # If there is a lower/equal-energy neighbor to merge with.
            if en < cg_ndata[nb]['energy']:
                continue
            barrier = cg_edata[(node, nb)]['saddle_energy'] - en
            if barrier < minh:
                break
        else: # it's a local minimum!
            cg_basin[node] = cg_basin.get(node, set())
            continue

        # it's a transition structure
        leqn = sorted([n for n in nbrs if cg_ndata[n]['energy'] <= en], 
                      key = lambda x: (cg_ndata[x]['energy'], x))

        # starting maximum barrier is just a tick lower than minh
        seen = set()
        reps, min_se = set(), en + minh - 1
        for nbr1 in leqn:
            assert ndata[nbr1]['energy'] <= en
            se1 = cg_edata[(node, nbr1)]['saddle_energy']
            if se1 < min_se: 
                # a neighbor with so far lowest saddle energy.
                reps, min_se = set([nbr1]), se1
            elif se1 == min_se: 
                # a neighbor with equally best saddle energy.
                reps.add(nbr1)
            for nbr2 in nbrs:
                if nbr1 == nbr2 or (nbr2, nbr1) in seen:
                    continue
                se2 = cg_edata[(nbr2, node)]['saddle_energy']
                se = max(se1, se2)
                if (nbr1, nbr2) in cg_edata:
                    se = min(se, cg_edata[(nbr1, nbr2)]['saddle_energy'])
                cg_edata[(nbr1, nbr2)] = {'saddle_energy': se}
                cg_edata[(nbr2, nbr1)] = {'saddle_energy': se}
                successors[nbr1].add(nbr2)
                successors[nbr2].add(nbr1)
                seen.add((nbr1, nbr2))

        # update the basins given the current representatives.
        if node in cg_basin:
            basin = cg_basin[node]
            del cg_basin[node]
        else:
            basin = set()
        for rep in reps: # Assume reps are lmins (for now).
            cg_basin[rep] = cg_basin.get(rep, set())
            cg_basin[rep].add(node)
            cg_basin[rep] |= basin

        # remove the node
        for nb in nbrs:
            del cg_edata[(nb, node)]
            del cg_edata[(node, nb)]
            successors[nb].remove(node)
        del cg_ndata[node]
        del successors[node]
    return cg_ndata, cg_edata, cg_basin

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Guide Landscpe construction                                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def get_guide_graph(seq, md, gnodes, gedges = None, tgn_folded = None):
    """
    """
    fc_empty = forbid_all_basepairs(seq, RNA.fold_compound(seq, md))
    tmp_gnodes = set(gnodes)
    assert len(tmp_gnodes) == len(gnodes)
    if tgn_folded is None:
        tgn_folded = {n: False for n in tmp_gnodes}
    else:
        assert all((n in tgn_folded) for n in gnodes)
    new_gnodes = []
    while True:
        # Do the constrained folding for all nodes.
        for con in gnodes:
            if not tgn_folded[con]:
                css, cen = mfe_constrained(seq, md, con)
                tgn_folded[con] = True
                tgn_folded[css] = True
                if css not in tmp_gnodes:
                    tmp_gnodes.add(css)
                    new_gnodes.append((css, cen))

        # Find all guide edges in the graph (only force input gedges).
        new_gedges = guiding_edge_search(tmp_gnodes, gedges)

        #NOTE: this may become too expensive. Adding all +1 edges?
        #print(f'{len(new_gedges)=}')
        guide_nbrs = {n: set() for n in tmp_gnodes}
        for (p, q) in new_gedges: 
            assert p != q
            guide_nbrs[p].add(q)
        for n in tmp_gnodes:
            for (p, q) in combinations(guide_nbrs[n], 2):
                assert (p, n) in new_gedges
                assert (n, q) in new_gedges
                dpq = get_bpd_cache(p, q)
                dpn = get_bpd_cache(p, n) 
                dnq = get_bpd_cache(n, q)
                if (p, q) in new_gedges:
                    # put an assert here
                    continue
                # NOTE: don't assert the line below,
                # the reason why two points are not connectec
                # can be a different node than current n.
                #assert max(dpn, dnq) <= dpq
                if dpn + dnq > dpq:
                    new_gedges.add((p, q))
                    new_gedges.add((q, p))
        #print(f'{len(new_gedges)=}')

        # Find guide nodes that should be added to the graph.
        add_gnodes = guiding_node_search(seq, md, tmp_gnodes, new_gedges, fc_empty)
        if not add_gnodes:
            break
        for (ss, en) in add_gnodes:
            assert ss not in tmp_gnodes
            assert ss not in tgn_folded
            tmp_gnodes.add(ss) 
            new_gnodes.append((ss, en))
            tgn_folded[ss] = False
    del fc_empty
    return new_gnodes, new_gedges

def guiding_edge_search(nodes, edges = None):
    """ Find all edges of the guiding neighborhood.

    Annoying case, cannot add pq at dpq = 8
     - p ..........((((.....((((.((.........)).)))).))))... 
     - q ...((((...)))).....((((.((.........)).))))........
     - i .(((......)))......((((.((.........)).))))........

    Args:
        edges (set, optional): Provide a set of edges. They will not be deleted!
    """
    if edges is None:
        edges = set()
    for p, q in combinations(nodes, 2):
        if (p, q) in edges:
            assert (q, p) in edges
            continue
        for i in nodes:
            if i == p or i == q:
                continue
            dpq = get_bpd_cache(p, q)
            dpi = get_bpd_cache(p, i)
            diq = get_bpd_cache(i, q)
            if max([dpi, diq]) < dpq:
            #if dpi + diq <= dpq:
                #print(f'Cannot add {p}, {q} at {dpq=} due to {i=} {dpi=} {diq=}.')
                break
            #dqp = get_bpd_i_cache(q, p)
            #dqi = get_bpd_i_cache(q, i)
            #dip = get_bpd_i_cache(i, p)
            #if max([dqi, dip]) < dqp:
            ##if dqi + dip <= dqp:
            #    #print(f'Cannot add {q}, {p} at {dqp=} due to {i=} {dqi=} {dip=}.')
            #    break

        else:
            edges.add((p, q))
            edges.add((q, p))
    return edges

def guiding_node_search(seq, md, nodes, edges, fc_empty, mind = 5):
    """ For every edge in the graph, find the MFE intersect.

    Args:
        ...
        fc_empty (fold compound): A fold compound where
            all base-pairs are forbidden.
        mind (int, optional): A minimum base-pair distance
            to do constrained folding. Defaults to 5.

    """
    seen = set()
    lmins = set()
    for (s1, s2) in edges:
        if get_bpd_cache(s1, s2) < mind:
            continue
        if (s2, s1) in seen:
            continue
        bps = get_basepairs([s1, s2])
        mss, mfe = mfe_intersect(seq, md, bps, fc_empty)
        if mss not in nodes:
            lmins.add((mss, mfe))
            css, cfe = mfe_constrained(seq, md, mss)
            if css not in nodes:
                lmins.add((css, cfe))
        seen.add((s1, s2))
    return lmins

def forbid_all_basepairs(seq, fc):
    for i in range(1, len(seq) + 1):
        for j in range(i + 4, len(seq) + 1):
            fc.hc_add_bp(i, j, RNA.CONSTRAINT_CONTEXT_NONE | RNA.CONSTRAINT_CONTEXT_NO_REMOVE)
    return fc

def get_basepairs(dotbrackets):
    bps = set()
    for db in dotbrackets:
        pt = RNA.ptable(db)
        bps |= set((i, j) for i, j in enumerate(pt[1:], 1) if j > i)
    return bps

def mfe_constrained(seq, md, dbcon):
    # TODO: is there a way to reset the fold compound?
    fc = RNA.fold_compound(seq, md)
    fc.constraints_add(dbcon, RNA.CONSTRAINT_DB_DEFAULT)# | RNA.CONSTRAINT_DB_ENFORCE_BP)
    mss, mfe = fc.mfe()
    return mss, int(round(mfe * 100))

def mfe_intersect(seq, md, bps, fc_empty = None):
    """ Return the MFE structure given allowed base pairs.

    Makes a temporary fold compound where all bases are forbidden and only the
    specified base-pairs are allowed. If fc is specified (for efficiency
    reasons only), then it is assumed that all base-pairs are already forbiden!
    The function will return the fold compund with forbidden base-pairs again.
    """
    fc_tmp = fc_empty is None 
    if fc_empty is None:
        fc_empty = forbid_all_basepairs(seq, RNA.fold_compound(seq, md))

    for bp in bps: # enable base pairs
        fc_empty.hc_add_bp(bp[0], bp[1], RNA.CONSTRAINT_CONTEXT_ALL_LOOPS | RNA.CONSTRAINT_CONTEXT_NO_REMOVE)
    mss, mfe = fc_empty.mfe()

    if fc_tmp:
        del fc_empty
    else:
        for bp in bps:
            fc_empty.hc_add_bp(bp[0], bp[1], RNA.CONSTRAINT_CONTEXT_NONE | RNA.CONSTRAINT_CONTEXT_NO_REMOVE)
    return mss, int(round(mfe * 100))
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Neighborhood processing                                                      #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
def local_flooding(fc, lmin, basinh = 2, rates = True, RT = None, k0 = None):
    """ Enumerate microstates within basin and the first steps out of the basin.

    Args:
        fc: The fold compound object.
        lmin: The local minimum representative secondary structure.
        basinh (kcal/mol): a local elevation higher than minh is considered a first step.
        rates (bool, optional): Get exit rates for all microtransitions.
    """
    macro = dict() # macro[ss] = [dG, P_ss] elevation and probability
    fstep = dict() # fstep[ss] = [(m1,k_ms), (m2,k_ms), ...] neighbors in macrostate and microrate

    dcbasinh = int(round(basinh*100))
    def local_nbrs(ss, refen):
        assert ss in macro # obvious, this was just assigned!
        for [nb, dcal] in get_neighbors(fc, ss):
            if dcal + refen >= dcbasinh:
                if nb in fstep:
                    fstep[nb].append((ss, round(dcal/100, 2)))
                else:
                    fstep[nb] = [(ss, round(dcal/100, 2))]
                continue
            if nb in macro:
                continue
            yield [nb, dcal]

    news = [(lmin, 0)]
    while news:
        newnews = set()
        for (newss, newen) in news:
            if newen < 0:
                #rlog.warning('Found a lower energy neighbor!')
                return None, newss
            assert newss not in fstep # this is going to be in macrostate!
            macro[newss] = [round(newen/100, 2), None] # elevation, probability
            for [nbr, dcal] in local_nbrs(newss, newen):
                newnews.add((nbr, newen+dcal))
        news = list(newnews)

    if rates:
        assert k0 is not None
        assert RT is not None
        # Alright, so we calculate the probabilities of states in the macrostate:
        x = np.array([m[0] for m in macro.values()])
        Z = sum(math.e**(-x/RT))
        for m in macro.keys():
            macro[m][1] = (math.e**(-macro[m][0]/RT))/Z 
        k_exit = 0
        for fs, ms in fstep.items():
            k_mj = 0
            P_mj = 0
            for [m, dG] in ms:
                k_ij = k0 * math.e**(-dG/RT)
                P_i = macro[m][1]
                k_mj += P_i * k_ij
                k_exit += P_i * k_ij
            fstep[fs] = [P_mj, k_mj]
        for fs, [P_mj, k_mj] in fstep.items():
            assert P_mj == 0
            fstep[fs][0] = fstep[fs][1]/k_exit
    return macro, fstep

def get_neighbors(fc, db = None, pt = None):
    """
    """
    assert (pt is None) or (db is None)
    if pt is None:
        pt = RNA.ptable(db)
    nbrs = []
    for move in fc.neighbors(pt):
        npt = list(pt)
        if move.is_removal():
            npt[-move.pos_3] = 0
            npt[-move.pos_5] = 0
            ndb = RNA.db_from_ptable(npt)
        elif move.is_insertion():
            npt[move.pos_3] = move.pos_5
            npt[move.pos_5] = move.pos_3
        else:
            rlog.warning(f"Are you using shift moves?")
            rlog.warning(f" shift = {move.is_shift()}")
            rlog.warning(f" pos3 = {move.pos_3}")
            rlog.warning(f" pos5 = {move.pos_5}")
            raise NotImplementedError('Are you using shift moves?')
        dG = fc.eval_move_pt(pt, move.pos_5, move.pos_3)
        if db:
            nbrs.append([RNA.db_from_ptable(npt), dG])
        else:
            nbrs.append([npt, dG])
    return nbrs

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Main stuff                                                                   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
def costruct(seq, cut = None):
    """ Translate between 'input' and 'internal' RNAcofold sequences """
    delim = '&'
    if delim in seq:
        cut = seq.index(delim) + 1
        table = str.maketrans(dict.fromkeys('&'))
        seq = seq.translate(table)
    elif cut and cut != -1:
        seq = seq[:cut-1] + delim + seq[cut-1:]
        cut = None
    return seq, cut

def parse_model_details(parser):
    """ ViennaRNA Model Details Argument Parser.  """
    model = parser.add_argument_group('ViennaRNA model details')

    model.add_argument("-T", "--temp", type = float, default = 37.0, 
        metavar = '<flt>',
        help = 'Rescale energy parameters to a temperature of temp C.')

    model.add_argument("-4", "--noTetra", action = "store_true",
        help = """Do not include special tabulated stabilizing 
        energies for tri-, tetra- and hexaloop hairpins.""")

    model.add_argument("-d", "--dangles", type = int, default = 2, 
        metavar = '<int>',
        help = """How to treat dangling end energies for bases adjacent to
        helices in free ends and multi-loops.""")

    model.add_argument("--noGU", action = "store_true",
        help = 'Do not allow GU/GT pairs.')

    model.add_argument("--noClosingGU", action = "store_true",
        help = 'Do not allow GU/GT pairs at the end of helices.')

    model.add_argument("-P", "--paramFile", action = "store", default = None,
        metavar = '<str>',
        help = """Read energy parameters from paramfile, instead of 
        using the default parameter set.""")

def main():
    import sys
    import argparse
    """ A wrapper for ViennaRNA findpath functions. """
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action='count', default = 0,
            help = "Verbose output, e.g. the folding pathway. (-vv increases verbosity.)")
    parser.add_argument("-w","--search-width", type = int, default = None, 
            help = "Adjust upper bound for findpath search.")
    parser.add_argument("-m","--max-energy", type = float, default = None,
            help = "Specify upper bound for barrier energy (kcal/mol).")
    parser.add_argument("--split", action = "store_true",
            help = "Split findpath into subproblems, if possible.")
    parser.add_argument("--minh", type = float, default = None,
            help = "Set a minimum barrier height for path flooding.")
    parse_model_details(parser)
    args = parser.parse_args()

    # TODO: Parse a subopt file.
    seq = None
    sss = dict()
    for e, line in enumerate(sys.stdin):
        if e == 0:
            seq = line.strip()
            continue
        l = line.split()
        if l[0]:
            sss[l[0]] = {'energy': int(round(float(l[1])*100))}
        else:
            print(f'Troubles with input line {e+1} {line=}')

    # Set model details.
    md = RNA.md()
    md.temperature = args.temp
    md.dangles = args.dangles
    md.special_hp = not args.noTetra
    md.logML = 0
    md.noGU = args.noGU
    md.noGUclosure = args.noClosingGU

    #import statprof
    #statprof.start()

    # Get the guide graph for all inputs:
    gnodes, gedges = get_guide_graph(seq, md, sss)

    if len(gnodes):
        print('Some important nodes seem to be missing in your input!')
        for gn in gnodes:
            print('', gn[0], gn[1])
            sss[gn[0]] = {'energy': gn[1]}

    print()
    for ge in gedges:
        print('', ge[0])
        print('', ge[1])
        print()

    ndata, edata = neighborhood_flooding((seq, md), sss, gedges, minh = args.minh)

    #statprof.stop()
    #statprof.display()

    #bpd = RNA.bp_distance(ss1, ss2)
    #rlog.info("Base-pair distance:", RNA.bp_distance(ss1, ss2))

    #if args.search_width is None:
    #    args.search_width = 2 * bpd
    #rlog.info(f"Search-width: {args.search_width}")
    #if args.verbose:
    #    print(f"   {seq}")
    #saddle = findpath_wrap(fc, ss1, ss2, 
    #        args.max_energy, 
    #        args.search_width, 
    #        cut = cut, 
    #        verbose = args.verbose)

    #if saddle is not None:
    #    e1 = round(fc.eval_structure(ss1), 2)
    #    barrier = saddle - e1
    #    print("Saddle: {:6.2f} kcal/mol | Barrier: {:6.2f} kcal/mol | Search width parameter: {}".format(saddle, barrier, args.search_width))
    #else:
    #    print("No valid saddle found below energy: {:6.2f} kcal/mol | Search width parameter: {}".format(args.max_energy, args.search_width))

    #if args.split:
    #    assert cut is None
    #    print('Findpath decomposition:')
    #    path, barrier = findpath_split(seq, ss1, ss2, md, th = 5, w = None)
    #    for x in path:
    #        print(x)

    #if args.minh is not None:
    #    assert args.split 
    #    assert cut is None
    #    ssmap = path_flooding(path, minh = int(args.minh*100))

    #    for si in sorted(ssmap):
    #        lm = ssmap[si]
    #        if isinstance(lm, list):
    #            assert len(lm) == 2
    #            print(*path[si], 'saddle')
    #        elif si == ssmap[si]:
    #            print(*path[si], 'lmin')
    #        else:
    #            print(*path[si])
 
if __name__ == '__main__':
    main()
    
