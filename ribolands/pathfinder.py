#!/usr/bin/env python
#
# ribolands.pathfinder
# 
# Utilities using the findpath algorithm.
#

import logging
rlog = logging.getLogger(__name__)

import sys
import RNA
import math
import argparse
import numpy as np
from itertools import islice, product, permutations

from ribolands.utils import make_loop_index, make_pair_table

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
    else :
        dist = RNA.bp_distance(s1, s2)
        BPD_CACHE[(s1, s2)] = dist
        return dist

def get_bp_change(seq, s1, s2):
    """ Extract (multiple) minimal regions where base-pairs change.

    Assume every common base-pair between s1 and s2 divides the direct path
    barrier estimation into two problems. This algorithm returns all pairs of
    substructures between s1 and s2 where base-pairs change, including the
    constant enclosing base-pairs (or adjacent dangling ends in the exterior
    loop case).

    TODO: This is a quick implementation, not particularly efficient.

    Args:
        seq (str): RNA sequence
        s1 (str): first structure
        s2 (str): second structure

    Returns:
        list[(seq, c1, c2),...]: A list of independent base-pair changes.
    """
    pt1 = make_pair_table(s1, base = 0, chars = list('.x'))
    pt2 = make_pair_table(s2, base = 0, chars = list('.x'))

    final_features = []

    feature_stack = []
    cse, cf1, cf2 = '', '', ''
    for e in range(0, len(pt1)):
        if pt1[e] == pt2[e] and pt1[e] != -1 : # a common base-pair
            if pt1[e] > e: 
                assert s1[e] == s2[e] and s1[e] == '('
                cse += seq[e]
                cf1 += s1[e]
                cf2 += s2[e]
                feature_stack.append([cse, cf1, cf2])
                cse, cf1, cf2 = seq[e], s1[e], s2[e]
            else:
                assert s1[e] == s2[e] and s1[e] == ')'
                cse += seq[e]
                cf1 += s1[e]
                cf2 += s2[e]
                if cf1 != cf2:
                    final_features.append((cse, cf1, cf2))
                [cse, cf1, cf2] = feature_stack.pop()
                cse += 'NNN' + seq[e]
                cf1 += 'xxx' + s1[e]
                cf2 += 'xxx' + s2[e]
        else:
            cse += seq[e]
            cf1 += s1[e]
            cf2 += s2[e]
    assert feature_stack == []

    # Exterior loop processing
    pf1 = make_pair_table(cf1, base = 0, chars = ['x','.'])
    pf2 = make_pair_table(cf2, base = 0, chars = ['x','.'])

    eloop = True
    record = 0
    xse, xf1, xf2 = '', '', ''
    pfiter = iter(range(0, len(pf1)))
    for e in pfiter:
        if pf1[e] == -1 and pf2[e] == -1 and eloop: # a single base-pair in the exerior loop
            assert cf1[e] == '.' and cf2[e] == '.'
            if xf1 != xf2:
                xse += cse[e]
                xf1 += cf1[e]
                xf2 += cf2[e]
                final_features.append((xse, xf1, xf2))
            xse, xf1, xf2 = cse[e], cf1[e], cf2[e]
            record = 0
        elif not eloop or pf1[e] != pf2[e]: # Start or keep recording.
            eloop = False
            record = max(record, pf1[e], pf2[e])
            xse += cse[e]
            xf1 += cf1[e]
            xf2 += cf2[e]
            if e == record:
                eloop = True
                record = 0
        elif pf1[e] == pf2[e]: # skip a common base-pair in the exterior loop
            xse += cse[e] + 'NNN' + cse[pf1[e]]
            xf1 += cf1[e] + 'xxx' + cf1[pf1[e]]
            xf2 += cf2[e] + 'xxx' + cf2[pf1[e]]
            next(islice(pfiter, pf1[e]-e-1, None))

    if xf1 != xf2:
        final_features.append((xse, xf1, xf2))

    return final_features

def apply_bp_change(seq, s1, s2, subseq, subs1, subs2):
    """Applies all bp changes of subs2 to s1.

    Note: s2 is only provided for sanity checks. s2 must be the output of this
        function, if subs2 is the only change from s1 to s2.

    Args:
        seq (str): Full length sequence
        s1 (str): Starting structure (containing subs1)
        s2 (str): Stop structure (containing subs2), can be None.
        subseq (str): subsequence of bpairs to be changed.
        subs1 (str): substructure to be changed from.
        subs2 (str): substructure to be changed to.

    Returns:
        [str]: The new structure.
    """
    subs = subseq.split('NNN')
    sb1s = subs1.split('xxx')
    sb2s = subs2.split('xxx')

    idxs = [[]] * len(subs)
    def find_indices(ssq, ss1, ss2):
        """Extract the indices where ssq + ss1 (+ ss2) match the parent structure.
        """
        slen = len(ssq)
        ixs = []

        idx = -1 
        while True:
            idx = seq.find(ssq, idx+1)
            if ss1 == s1[idx:idx+slen] and (s2 is None or ss2 == s2[idx:idx+slen]):
                ixs.append(idx)
            if idx == -1:
                break
        return ixs 

    for e in range(len(subs)):
        idxs[e] = find_indices(subs[e], sb1s[e], sb2s[e])

    def get_news_iter():
        """ Iterate over all new structures that can be composed from the indices.
        """
        li = make_loop_index(s1)
        for x in product(*idxs):
            # There is some excess work done here, 
            # because not all products are valid combiations ...
            loops = []
            news = s1
            assert len(x) == len(subs)
            idx = x[0]
            alen = len(subs[0])
            news = news[0:idx] + sb2s[0] + news[idx+alen:]
            for a, b in enumerate(range(1, len(x))):
                alen = len(subs[a])
                blen = len(subs[b])
                aid = x[a]
                bid = x[b]
                if aid + alen >= bid:
                    break
                elif li[aid] != li[bid+blen-1] and li[bid] != li[aid+alen-1]:
                    # At least one of them has to match ...
                    break
                else:
                    if li[aid] == li[bid+blen-1]:
                        loops.append(li[aid])
                    if li[bid] == li[aid+alen-1]:
                        loops.append(li[bid])
                    news = news[0:bid] + sb2s[b] + news[bid+blen:]
            else:
                yield news, min(loops)
   
    if not all(len(x)==1 for x in idxs):
        # do the more expensive stuff
        results = set(get_news_iter())
        if len(results) > 1:
            # NOTE: This is full of hacks and special cases ... 
            # very, very bad style!
            exteriors = []
            for new, ml in results:
                if ml == 0:
                    exteriors.append(new)
            if len(exteriors) > 1:
                rlog.warning('Ambiguous base-pair change!')
                rlog.warning(f"SQ: {seq}, {subseq}")
                rlog.warning(f"S1: {s1}, {subs1}")
                rlog.warning(f"S2: {s2}, {subs2}")
                rlog.warning(f"Indices: {idxs}")
                for new, ml in results:
                    rlog.warning(f" - {new}, {ml}")
                raise PathfinderError('Ambiguous base-pair change!')
        news = list(results)[0][0] # list index out of range?
    else:
        news = s1
        for e, site in enumerate(idxs):
            idx = site[0]
            slen = len(subs[e])
            news = news[0:idx] + sb2s[e] + news[idx+slen:]
    return news

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Cache to look-up findpath barriers                                           #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# FPATH_CACHE[(seq, ss1, ss2)] = (round(sE - e1, 2), round(e2 - e1, 2), maxbar, fpw)
# FPATH_CACHE[(seq, ss1, ss2)] = (None, [(cs1, ss1, a), (cs1, ss1, ss3), ...])
FPATH_CACHE = {}

def clear_fpath_cache():
    global FPATH_CACHE
    FPATH_CACHE = {}

def get_fpath_cache(seq, ss1, ss2, md, maxbar = None, fpw = None):
    """ Store/update/return the cache for a given findpath calculation.

    Note: This routine accepts any input for ss1 and ss2. There is no restriction
        on whether it is a minimal change or not. However, as a fold compound is
        generated every time this function is called, it is more suitable for a 
        context where the input sequence changes, such as when regions of minimal
        change are given.

    TODO: Unify the usage in combination with get_fpath_flooding_cache().

    Args:
        seq (str): Nucleic acid sequence.
        ss1 (str): Start secondary structure.
        ss2 (str): Stop secondary structure.
        md (obj): ViennaRNA model details object.
        maxbar (float, optional): The maximum acceptable barrier. If maxbar is None or
            float('inf'), then a barrier must be found. Defaults to None.
        fpw (int, optional): The findpath width parameter. If not specified, it
            is calculated as twice the base-pair distance.  Defaults to None.

    Returns:
        FPATH_CACHE[(seq, ss1, ss2)]
    """
    global FPATH_CACHE

    if maxbar is None:
        maxbar = float('inf')
    if fpw is None:
        fpw = 2 * get_bpd_cache(ss1, ss2) 

    def fp_wrap(seq, ss1, ss2, md, maxbar, fpw):
        # Call findpath and store metadata.
        fc = RNA.fold_compound(seq, md)
        e1 = round(fc.eval_structure(ss1), 2)
        e2 = round(fc.eval_structure(ss2), 2)

        # Default mode: do findpath and store the results.
        if maxbar == float('inf'):
            dcal_sE = fc.path_findpath_saddle(ss1, ss2, width = fpw)
        else:
            dcal_bound = int(round((maxbar + e1) * 100))
            dcal_sE = fc.path_findpath_saddle(ss1, ss2, maxE = dcal_bound, width = fpw)

        sE = float(dcal_sE)/100 if dcal_sE is not None else float('inf')
        FPATH_CACHE[(seq, ss1, ss2)] = (round(sE - e1, 2), round(e2-e1, 2), maxbar, fpw)
        FPATH_CACHE[(seq, ss2, ss1)] = (round(sE - e2, 2), round(e1-e2, 2), maxbar, fpw)
        del fc
 
    if (seq, ss1, ss2) not in FPATH_CACHE:
        fp_wrap(seq, ss1, ss2, md, maxbar, fpw)
    else:
        (obar, ogain, omax, ofpw) = FPATH_CACHE[(seq, ss1, ss2)]
        if obar == float('inf') and (maxbar > omax or fpw > ofpw): # retry with a higher bound.
            fp_wrap(seq, ss1, ss2, md, maxbar, fpw)
        elif obar != float('inf') and (obar + 0.01 <= maxbar and fpw > ofpw):
            fp_wrap(seq, ss1, ss2, md, maxbar, fpw)
        assert obar >= FPATH_CACHE[(seq, ss1, ss2)][0]
        assert ogain == FPATH_CACHE[(seq, ss1, ss2)][1]
    return FPATH_CACHE[(seq, ss1, ss2)]

def path_flooding(path, minh):
    """Use flooding algorithm to determine local minima on a folding path.
    
    Identifies the lowest energy of the local minimum, one representative
    stucture and "other structures" associated with that minimum. Beware that
    the "other structures" can contain saddle components when the function
    is applied to paths with degenerate saddle components.

    Note: Local minima on a folding path are only local minima with respect to
    that path, not with respect to the full neighborhood.

    Args:
        path (list): A list of tuples (ss, en).
        minh (float): minimal height of an energy barrier separatig two basins.

    Returns:
        [dict, dict]: properties of the local minima.
    """
    indexpath = [[ss, en, step] for step, [ss, en] in enumerate(path)]
    sconf = sorted(indexpath, key = lambda x: x[1]) # by energy

    def get_ss(step):
        return path[step][0]

    lmins = dict() # lmins[ss] = [en, set(steps)]
    ssmap = dict() # ssmap[step] = step

    for e, (ss, en, step) in enumerate(sconf):
        if step == 0 or step == len(path)-1: # edge cases
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
                if en - highen < minh - 0.0001: # merge higher minimum into lower one.
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
    return lmins, ssmap

def get_fpath_flooding_cache(seq, ss1, ss2, md, minh, maxbar = None, fpwf = 2):
    """ Chop a direct folding path into multiple (independent?) steps.

    Note: This routine uses the same global FPATH_CACHE as get_fpath_cache, but 
        with a different format. They must not be used together!

    Args:
        seq (str): Nucleic acid sequence.
        ss1 (str): Start secondary structure.
        ss2 (str): Stop secondary structure.
        md (obj): ViennaRNA model details object.
        maxbar (float, optional): The maximum acceptable barrier. If maxbar is None or
            float('inf'), then a barrier will be found.  Defaults to None.
        fpwf (int, optional): The findpath width factor. For every findpath call,
            the width is calculated as fpwf * base-pair distance.  Defaults to 2.

    Returns:
        FPATH_CACHE[(seq, ss1, ss2)]
    """

    global FPATH_CACHE

    if maxbar is None:
        maxbar = float('inf')

    def findpath_get_path(seq, ss1, ss2, md, maxbar):
        # Call findpath and store metadata.
        fpw = fpwf * get_bpd_cache(ss1, ss2) 
        fc = RNA.fold_compound(seq, md)
        e1 = round(fc.eval_structure(ss1), 2)
        e2 = round(fc.eval_structure(ss2), 2)
        if maxbar == float('inf'):
            path = fc.path_findpath(ss1, ss2, width = fpw)
        else:
            dcal_bound = int(round((maxbar + e1) * 100))
            path = fc.path_findpath(ss1, ss2, maxE = dcal_bound, width = fpw)
        del fc # ?
        return [(step.s, step.en) for step in path] if len(path) else None

    def decompose_path(seq, ss1, ss2, firsttime = False):
        # Find a set of orthogonal pathways.
        # NOTE: We always start and stop with decomposition.
        FPATH_CACHE[(seq, ss1, ss2)] = [None, []]
        # Assume we will have to refer to some subproblem ...
        chopped = get_bp_change(seq, ss1, ss2)
        rlog.debug(f'My chopped path: {chopped}')
        if firsttime or len(chopped) > 1 or chopped[0][0] != seq:
            # Always do at least one cycle of flooding.
            # Iterate over a set of orthogonal of steps.
            FPATH_CACHE[(seq, ss1, ss2)][1] = set(chopped)
            for (msq, ms1, ms2) in chopped:
                assert len(ms1) == len(ms2)
                lminp = flood_path(msq, ms1, ms2)
                # iterate over a sequence of steps.
                for (ls1, ls2) in lminp:
                    assert len(ls1) == len(ls2)
                    decompose_path(msq, ls1, ls2)
        else: # stop condition
            assert len(chopped) == 1
            (msq, ms1, ms2) = chopped[0]
            assert [msq, ms1, ms2] == [seq, ss1, ss2]
            # Do the findpath calculation.
            if FPATH_CACHE[(msq, ms1, ms2)][0] is None:
                del FPATH_CACHE[(msq, ms1, ms2)] # delete the placeholder key
            _ = get_fpath_cache(msq, ms1, ms2, md, maxbar) # caching / final entry

    def flood_path(msq, ms1, ms2):
        # Return a sequence of start/end structures along all path minima (using minh).
        if (msq, ms1, ms2) in FPATH_CACHE and FPATH_CACHE[(msq, ms1, ms2)][0] is not None:
            pminpath = [(ms1, ms2)]
        elif (msq, ms1, ms2) in FPATH_CACHE and isinstance(FPATH_CACHE[(msq, ms1, ms2)][1], list): 
            pminpath = FPATH_CACHE[(msq, ms1, ms2)][1]
        else:
            path = findpath_get_path(msq, ms1, ms2, md, maxbar)
            if path is None: 
                raise NotImplementedError
            lmins, ssmap = path_flooding(path, minh)
            
            pminpath = []
            rpminpath = []
            ss1 = path[0][0]
            for ss, en in path[1:-1]: 
                if ss in lmins:
                    pminpath.append((ss1, ss))
                    rpminpath.insert(0, (ss, ss1))
                    ss1 = ss
            pminpath.append((ss1, path[-1][0]))
            rpminpath.insert(0, (path[-1][0], ss1))
            FPATH_CACHE[(msq, ms1, ms2)] = (None, pminpath)
            FPATH_CACHE[(msq, ms2, ms1)] = (None, rpminpath)
        return pminpath

    if (seq, ss1, ss2) not in FPATH_CACHE:
        decompose_path(seq, ss1, ss2, firsttime = True)
    elif FPATH_CACHE[(seq, ss1, ss2)][0] is None:
        pass
    else:
        (obar, ogain, omax, ofpw) = FPATH_CACHE[(seq, ss1, ss2)]
        fpw = fpwf * get_bpd_cache(ss1, ss2) 
        if obar == float('inf') and (maxbar > omax or fpw > ofpw): # retry with a higher bound.
            decompose_path(seq, ss1, ss2, firsttime = True)
        elif obar != float('inf') and (obar + 0.01 <= maxbar and fpw > ofpw):
            decompose_path(seq, ss1, ss2, firsttime = True)
        assert obar >= FPATH_CACHE[(seq, ss1, ss2)][0]
        assert ogain == FPATH_CACHE[(seq, ss1, ss2)][1]
    return FPATH_CACHE[(seq, ss1, ss2)]

def show_flooded_prime_path(seq, ss1, ss2):
    """ Returns ONE prime path, but there may be many.
    The cache must be filled before!
    """
    ppmlist = []
    def my_ppm(seq, ss1, ss2, startseq, start, stop):
        lmp = get_fpath_flooding_cache(seq, ss1, ss2, None, None)
        if lmp[0] is None:
            if isinstance(lmp[1], set):
                backup = start
                for perm in permutations(lmp[1]):
                    start = backup
                    for (lsq, pm1, pm2) in perm: 
                        my_ppm(lsq, pm1, pm2, startseq, start, stop)
                        start = apply_bp_change(startseq, start, stop, lsq, pm1, pm2)
            else:
                for (pm1, pm2) in lmp[1]: 
                    my_ppm(seq, pm1, pm2, startseq, start, None)
                    start = apply_bp_change(startseq, start, None, seq, pm1, pm2)
        else:
            stop = apply_bp_change(startseq, start, stop, seq, ss1, ss2)
            revp = get_fpath_flooding_cache(seq, ss2, ss1, None, None)
            ppmlist.append([start, stop, lmp[0], revp[0]])
            return start, stop, lmp[0], revp[0]
    my_ppm(seq, ss1, ss2, seq, ss1, ss2)
    return ppmlist

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Neighborhood processing                                                      #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def local_flooding(fc, lmin, basinh = 2, rates = True, RT = None, k0 = None):
    """ Enumerate microstates within basin and the first steps out of the basin.

    Args:
        fc: The fold compound object.
        lmin: The local minimum representative secondary structure.
        basinh: a local elevation higher than minh is considered a first step.
        rates (bool, optional): Get exit rates for all microtransitions.
    """
    macro = dict() # macro[ss] = [dG, P_ss] elevation and probability
    fstep = dict() # fstep[ss] = [(m1,k_ms), (m2,k_ms), ...] neighbors in macrostate and microrate

    def local_nbrs(ss, refen):
        assert ss in macro # obvious, this was just assigned!
        for [nb, dG] in get_neighbors(fc, ss):
            dG = dG/100
            if dG + refen >= basinh + 0.00001:
                if nb in fstep:
                    fstep[nb].append((ss, dG))
                else:
                    fstep[nb] = [(ss, dG)]
                continue
            if nb in macro:
                continue
            yield [nb, dG]

    news = [(lmin, 0)]
    while news:
        newnews = set()
        for (newss, newen) in news:
            assert newen >= 0 # elevation cannot be < 0
            macro[newss] = [newen, None] # elevation, probability
            if newss in fstep:
                print(newss, newen)
            assert newss not in fstep # this is going to be in macrostate!
            for [nbr, dG] in local_nbrs(newss, newen):
                newnews.add((nbr, newen+dG))
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
        p_tot = 0
        for fs, [P_mj, k_mj] in fstep.items():
            assert P_mj == 0
            fstep[fs][0] = fstep[fs][1]/k_exit
    return macro, fstep

def get_neighbors(fc, db = None, pt = None):
    """
    """
    if pt is None:
        pt = RNA.ptable(db)
    else:
        assert db == None
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
# Findpath                                                                     #
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

def findpath_wrap(fc, s1, s2, maxE, fpath, cut = None, verbose = False):
    """ A wrapper for ViennaRNA findpath functions. """
    sE = None
    if verbose:
        if maxE:
            dcal_bound = int(round(maxE * 100))
            path = fc.path_findpath(s1, s2, 
                    maxE = dcal_bound, width = fpath)
        else:
            path = fc.path_findpath(s1, s2, width = fpath)

        if len(path):
            sE = round(fc.eval_structure(s1), 2)
            for e, step in enumerate(path):
                ss, _ = costruct(step.s, cut)
                print("{:2d} {:s} {:6.2f}".format(e, ss, step.en))
                if step.en > sE: sE = step.en
    else :
        if maxE:
            dcal_bound = int(round(maxE * 100))
            dcal_sE = fc.path_findpath_saddle(s1, s2, 
                    maxE = dcal_bound, width = fpath)
        else :
            dcal_sE = fc.path_findpath_saddle(s1, s2, width = fpath)
        if dcal_sE is not None:
            sE = float(dcal_sE)/100
    return sE


def parse_model_details(parser):
    """ ViennaRNA Model Details Argument Parser.  """
    model = parser.add_argument_group('ViennaRNA model details')

    model.add_argument("-T", "--temp", type = float, default = 37.0, metavar = '<flt>',
        help = 'Rescale energy parameters to a temperature of temp C.')

    model.add_argument("-4", "--noTetra", action = "store_true",
        help = 'Do not include special tabulated stabilizing energies for tri-, tetra- and hexaloop hairpins.')

    model.add_argument("-d", "--dangles", type = int, default = 2, metavar = '<int>',
        help = 'How to treat "dangling end" energies for bases adjacent to helices in free ends and multi-loops.')

    model.add_argument("--noGU", action = "store_true",
        help = 'Do not allow GU/GT pairs.')

    model.add_argument("--noClosingGU", action = "store_true",
        help = 'Do not allow GU/GT pairs at the end of helices.')

    model.add_argument("-P", "--paramFile", action = "store", default = None, metavar = '<str>',
        help = 'Read energy parameters from paramfile, instead of using the default parameter set.')

def main():
    """ A wrapper for ViennaRNA findpath functions. """
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action='count', default = 0,
            help = "Verbose output, e.g. the folding pathway. (-vv increases verbosity.)")
    parser.add_argument("-w","--search-width", type = int, default = None, 
            help = "Adjust upper bound for findpath search.")
    parser.add_argument("--minh", type = float, default = None,
            help = "Set a minimum barrier height for path flooding.")
    parser.add_argument("-m","--max-energy", type = float, default = None,
            help = "Specify upper bound for barrier energy (kcal/mol).")
    parser.add_argument("--primepaths", action = "store_true",
            help = "Decompose findpath into primepaths.")
    parse_model_details(parser)
    args = parser.parse_args()

    # Verbose level 0: show output
    # Verbose level 1: show path
    # Verbose level 2: show rlog.info
    # Verbose level 3: show rlog.debug

    for e, line in enumerate(sys.stdin):
        if e == 0 :
            seq = line.strip()
        elif e == 1 :
            ss1, cut = costruct(line.strip())
        elif e == 2 :
            ss2, cut_ = costruct(line.strip())

    if cut != cut_:
        raise PathfinderError('Inconsistent cut-points.')

    # Set model details.
    md = RNA.md()
    md.temperature = args.temp
    md.dangles = args.dangles
    md.logML = 0
    md.special_hp = not args.noTetra
    md.noGU = args.noGU
    md.noGUclosure = args.noClosingGU
    fc = RNA.fold_compound(seq, md)

    bpd = RNA.bp_distance(ss1, ss2)
    rlog.info("Base-pair distance:", RNA.bp_distance(ss1, ss2))

    if args.search_width is None:
        args.search_width = 2 * bpd
    rlog.info(f"Search-width: {args.search_width}")
    if args.verbose:
        print(f"   {seq}")
    saddle = findpath_wrap(fc, ss1, ss2, 
            args.max_energy, 
            args.search_width, 
            cut = cut, 
            verbose = args.verbose)

    if saddle is not None:
        e1 = round(fc.eval_structure(ss1), 2)
        barrier = saddle - e1
        print("Saddle: {:6.2f} kcal/mol | Barrier: {:6.2f} kcal/mol | Search width parameter: {}".format(saddle, barrier, args.search_width))
    else:
        print("No valid saddle found below energy: {:6.2f} kcal/mol | Search width parameter: {}".format(args.max_energy, args.search_width))

    if args.primepaths:
        print('Findpath decomposition:')
        features = get_bp_change(seq, ss1, ss2)
        for (cseq, cs1, cs2) in features:
            print(cseq, cs1, cs2, get_fpath_cache(cseq, cs1, cs2, md))
        clear_fpath_cache()

    if args.minh is not None:
        rlog.info(f'Cache -- minh = {args.minh}:')
        get_fpath_flooding_cache(seq, ss1, ss2, md, args.minh, maxbar = None, fpwf = 2)
        pp = show_flooded_prime_path(seq, ss1, ss2)
        for [a, b, fw, rv] in pp:
            print(f"{a} -> {b} [fw = {fw}, rv = {rv}]")

if __name__ == '__main__':
    main()

