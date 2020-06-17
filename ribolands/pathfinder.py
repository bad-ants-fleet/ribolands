import RNA
import math
import numpy as np
import ribolands as ril


""" Utils around the findpath algorithm:
"""

def path_flooding(path, minh = 1):
    """Use flooding algorithm to determine local minima in a folding path.
    
    Identifies the lowest energy of the local minimum, one representative
    stucture and "other structures" associated with that minimum. Beware that
    the "other structures" can contain saddle components when the function
    is applied to paths with degenerate saddle components.

    Args:
        path (list): A list of tuples (ss, en).

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
                if left == right: # boring, merge!
                    lms, lme = path[left]
                    assert en >= lme
                    lmins[lms][1].add(step)
                    ssmap[step] = left
                else: # yay, saddle!
                     lmsl, lmel = path[left]
                     lmsr, lmer = path[right]
                     [lower, higher] = [left, right] if lmel < lmer else [right, left]
                     [lowlm, highlm] = [lmsl, lmsr] if lmel < lmer else [lmsr, lmsl]
                     [lowen, highen] = [lmel, lmer] if lmel < lmer else [lmer, lmel]
                     if en - highen <= minh: # merge higher minimum into lower one.
                         lmins[lowlm][1] |= lmins[highlm][1] | set([step]) 
                         for hlm in lmins[highlm][1]:
                             ssmap[hlm] = lower
                         ssmap[step] = lower
                         del lmins[highlm]
                     else: # keep them separated, but what do we do with ss?
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

FPATH_CACHE = {}

def get_fpath_cache(seq, cs1, cs2, md, maxbar, fpw):
    """Store and update the barrier for a minimal region enclosed by static base-pairs.

    TODO: develop an indirect-barrier concept: Under which circumstances can
        you find a better indirect path for a single orthogonal region?

        Approach: Since we know sequence, start structure and end structure, we
        could ask if we have previously found a different path cs1 -> ts ->
        cs2, but then we would have to adapt the datastructure.

    Args:
        maxbar (float): The maximum acceptable barrier. If maxbar is None or
            float('inf'), then a barrier will be found.

    """
    global FPATH_CACHE

    if maxbar is None:
        maxbar = float('inf')
    assert fpw >= get_bpd_cache(cs1, cs2) 

    def fp_wrap(seq, cs1, cs2, md, maxbar, fpw):
        # Call findpath and store metadata.
        fc = RNA.fold_compound(seq, md)
        e1 = round(fc.eval_structure(cs1), 2)
        e2 = round(fc.eval_structure(cs2), 2)

        if maxbar == float('inf'):
            dcal_sE = fc.path_findpath_saddle(cs1, cs2, width = fpw)
        else:
            dcal_bound = int(round((maxbar + e1) * 100))
            dcal_sE = fc.path_findpath_saddle(cs1, cs2, maxE = dcal_bound, width = fpw)

        sE = float(dcal_sE)/100 if dcal_sE is not None else float('inf')

        FPATH_CACHE[(seq, cs1, cs2)] = (round(sE - e1, 2), round(e2-e1, 2), maxbar, fpw)
        FPATH_CACHE[(seq, cs2, cs1)] = (round(sE - e2, 2), round(e1-e2, 2), maxbar, fpw)
        del fc
 
    if (seq, cs1, cs2) not in FPATH_CACHE:
        fp_wrap(seq, cs1, cs2, md, maxbar, fpw)
    else:
        (obar, ogain, omax, ofpw) = FPATH_CACHE[(seq, cs1, cs2)]
        if obar == float('inf') and (maxbar > omax or fpw > ofpw): # retry with a higher bound.
            fp_wrap(seq, cs1, cs2, md, maxbar, fpw)
        elif obar != float('inf') and (obar + 0.01 <= maxbar and fpw > ofpw):
            fp_wrap(seq, cs1, cs2, md, maxbar, fpw)

        assert obar >= FPATH_CACHE[(seq, cs1, cs2)][0]
        assert ogain == FPATH_CACHE[(seq, cs1, cs2)][1]

    return FPATH_CACHE[(seq, cs1, cs2)]

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
    pt1 = ril.make_pair_table(s1, base = 0)
    pt2 = ril.make_pair_table(s2, base = 0)

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
    pf1 = ril.make_pair_table(cf1, base = 0, chars = ['x','.'])
    pf2 = ril.make_pair_table(cf2, base = 0, chars = ['x','.'])

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


def local_flooding(fc, lmin, minh = 2, rates = True, RT = None, k0 = None):
    """
    """
    macro = dict() # macro[ss] = [dG, P_ss] elevation and probability
    fstep = dict() # fstep[ss] = [(m1,k_ms), (m2,k_ms), ...] neighbors in macrostate and microrate

    def local_nbrs(ss, refen):
        assert ss in macro
        for [nb, dG] in get_neighbors(fc, ss):
            dG /= 100
            if dG + refen > minh:
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
            macro[newss] = [newen, None]
            assert newss not in fstep
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
            print(move.is_shift())
            print(move.pos_3)
            print(move.pos_5)
            raise NotImplementedError('shift moves?')

        dG = fc.eval_move_pt(pt, move.pos_5, move.pos_3)
        if db:
            nbrs.append([RNA.db_from_ptable(npt), dG])
        else:
            nbrs.append([npt, dG])

    return nbrs

