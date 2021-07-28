#!/usr/bin/env python3

import RNA
import unittest
from itertools import combinations

from ribolands.utils import make_pair_table
from ribolands.pathfinder import (get_bpd_cache, 
                                  get_bpd_i_cache,
                                  findpath_split,
                                  common_basepairs,
                                  common_exterior_bases, 
                                  split_struct,
                                  merge_struct,
                                  local_flooding,
                                  path_flooding, 
                                  edge_flooding,
                                  init_findpath_max,
                                  findpath_max,
                                  get_guide_graph,
                                  guiding_edge_search,
                                  guiding_node_search,
                                  forbid_all_basepairs,
                                  get_basepairs,
                                  mfe_intersect,
                                  neighborhood_flooding,
                                  top_down_coarse_graining)

from ribolands.parser import parse_barriers

SKIP = False

@unittest.skipIf(SKIP, "skipping tests")
class MaxPathTests(unittest.TestCase):
    #def test_cache(self):
    #    search_width_multiplier = 4
    #    mp = True
    #     
    #    sequence = 'UCCGACAUUAAGACACACCAGGGUCUCGUAUCCCUAGGGUAAGGUACGCGCGGACCGGCCAAUCGGGUAUUGCUGCAAACUAUGGCAAUAGUGCAUAGGUUCAGACGAAGUACGGGUGGAUAUUUGUAGCCAGUAUGCUGGGUCUCCGGG'
    #    fp = maxpath.findpath_class(sequence, mp)
    #    
    #    s1       = '((((..........((.((((((........)))).))))..........))))((((...(((.(((((.((((((((((((.(((....)))))))(((((..((.....))..))))).))))))))....))))).)))..)))).'
    #    s2       = '((((..........((.((((((........)))).))))..........))))((((....((((((((((((((((((((((((....)).)))))(((((..((.....))..))))).)))))))).))))).))))....)))).'
    #    result = fp.init(s1, s2, search_width_multiplier)
    #    print(result)
    #    
    #    s1       = '((((..........((.((((((........)))).))))..........))))((((...(((.(((((.((((((((((((.(((....)))))))(((((..((.....))..))))).))))))))....))))).)))..)))).'
    #    s2       = '((((....((....((.((((((........)))).))))....))....))))((((....((((((((((((((((((((((((....)).)))))(((((..((.....))..))))).))))))).)))))).))))....)))).'
    #    result = fp.init(s1, s2, search_width_multiplier)
    #    print(result)


    def test_findpath_split_02(self):
        seq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCG"
        ss1 = "..((....(((........(((.....)))....)))..(((........)))..))." 
        ss2 = ".((((....))))..(((.(.((...)).)..)))......(((.......)))...."
        #print(f'\n{seq}\n{ss1}\n{ss2}')
        path, barrier = findpath_split(seq, ss1, ss2, RNA.md(), th = 1)
        assert barrier == 200 # may change with a new findpath version?

    def test_mfe_intersect(self):
        seq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCG"
        ss1 = "..((....(((........(((.....)))....)))..(((........)))..))." 
        ss2 = ".((((....))))..(((.(.((...)).)..)))......(((.......)))...."
        bps = get_basepairs([ss1, ss2])
        mss, mfe = mfe_intersect(seq, RNA.md(), bps)
        assert mss != ss1
        assert mss != ss2

    def test_mfe_intersect(self):
        seq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCG"
        ss1 = "..((....(((........(((.....)))....)))..(((........)))..))." 
        ss2 = ".((((....))))..(((.((((...))))..)))....(((........)))....."
        bps = get_basepairs([ss1, ss2])
        mss, mfe = mfe_intersect(seq, RNA.md(), bps)
        assert mss == ss2
 
    def test_mfe_intersect(self):
        seq = 'UCCGACAUUAAGACACACCAGGGUCUCGUAUCCCUAGGGUAAGGUACGCGCGGACCGGCCAAUCGGGUAUUGCUGCAAACUAUGGCAAUAGUGCAUAGGUUCAGACGAAGUACGGGUGGAUAUUUGUAGCCAGUAUGCUGGGUCUCCGGG'
        ss1 = '((((..........((.((((((........)))).))))..........))))((((...(((.(((((.((((((((((((.(((....)))))))(((((..((.....))..))))).))))))))....))))).)))..)))).'
        ss2 = '((((..........((.((((((........)))).))))..........))))((((....((((((((((((((((((((((((....)).)))))(((((..((.....))..))))).)))))))).))))).))))....)))).'
        ss3 = '((((..........((.((((((........)))).))))..........))))((((....((((((((((((((((((((((.((....)))))))(((((..((.....))..))))).)))))))).))))).))))....)))).'
        bps = get_basepairs([ss1, ss2, ss3])
        mss, mfe = mfe_intersect(seq, RNA.md(), bps)
        assert mss == ss3

@unittest.skipIf(SKIP, "skipping tests")
class FindpathTests(unittest.TestCase):
    def test_split_merge_01(self):
        seq = 'GCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUU'
        ss1 = '(......((((((((......)))).......))))..)'
        ss2 = '((((.......((((......)))).......)))...)'
        pt1 = make_pair_table(ss1, base = 0)
        pt2 = make_pair_table(ss2, base = 0)

        ceb = []
        assert list(common_exterior_bases(pt1, pt2)) == ceb
        for x in common_exterior_bases(pt1, pt2):
            assert seq == merge_struct(*split_struct(seq, x, None), x, None)
            assert ss1 == merge_struct(*split_struct(ss1, x, None), x, None)
            assert ss2 == merge_struct(*split_struct(ss2, x, None), x, None)

        cbp = [(0, 38), (11, 24), (12, 23), (13, 22), (14, 21)]
        assert list(common_basepairs(pt1, pt2)) == cbp
        for x in common_basepairs(pt1, pt2):
            assert seq == merge_struct(*split_struct(seq, *x), *x)
            assert ss1 == merge_struct(*split_struct(ss1, *x), *x)
            assert ss2 == merge_struct(*split_struct(ss2, *x), *x)

    def test_split_merge_02(self):
        seq = 'AAAGCCGCCUUAAACGGGUAUUGGUACCNNNGGCAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUAC'
        ss1 = '...((((.(((....)))...)))).((...))...((((.(((..((.(((((......)))))..))..)))...))))..'
        ss2 = '....((((((.....))))...))..((...))....((((((.......((((......)))).......))))...))...'
        pt1 = make_pair_table(ss1, base = 0)
        pt2 = make_pair_table(ss2, base = 0)

        ceb = [0, 1, 2, 25, 33, 34, 35, 81, 82]
        assert list(common_exterior_bases(pt1, pt2)) == ceb
        for x in common_exterior_bases(pt1, pt2):
            assert seq == merge_struct(*split_struct(seq, x, None), x, None)
            assert ss1 == merge_struct(*split_struct(ss1, x, None), x, None)
            assert ss2 == merge_struct(*split_struct(ss2, x, None), x, None)

        cbp = [(4, 23), (5, 22), (26, 32), (27, 31), (37, 79), (38, 78), (50, 63), (51, 62), (52, 61), (53, 60)]
        assert list(common_basepairs(pt1, pt2)) == cbp
        for x in common_basepairs(pt1, pt2):
            assert seq == merge_struct(*split_struct(seq, *x), *x)
            assert ss1 == merge_struct(*split_struct(ss1, *x), *x)
            assert ss2 == merge_struct(*split_struct(ss2, *x), *x)

    def test_findpath_split_00(self):
        seq = 'GCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUU'
        ss1 = '(......((((((((......)))).......))))..)'
        ss2 = '((((.......((((......)))).......)))...)'
        # Set model details.
        path, barrier = findpath_split(seq, ss1, ss2, RNA.md())
        assert path == [('(......((((((((......)))).......))))..)', 400),
                        ('(.......(((((((......)))).......)))...)', 600),
                        ('(.......((.((((......))))........))...)', 860),
                        ('(........(.((((......))))........)....)', 1210),
                        ('(..........((((......)))).............)', 680),
                        ('((.........((((......)))).........)...)', 770),
                        ('(((........((((......))))........))...)', 430),
                        ('((((.......((((......)))).......)))...)', 280)]
        assert barrier == max((en for ss, en in path)) - path[0][1]

    def test_findpath_split_01(self):
        seq = 'AAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUACAAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUAC'
        ss1 = '...((((.(((..((.(((((......)))))..))..)))...))))...............((((.(((..((.(((((......)))))..))..)))...))))............'
        ss2 = '...(((((((.......((((......)))).......))))...)))...............(((((((.......((((......)))).......))))...)))............'
       #ss2 = '....((((((.......((((......)))).......))))...)).................((((((.......((((......)))).......))))...)).............'
        vrna_md = RNA.md()
        #print(f'\n{seq}\n{ss1}\n{ss2}')
        path, barrier = findpath_split(seq, ss1, ss2, vrna_md, th = 99)
        assert path == [('...((((.(((..((.(((((......)))))..))..)))...))))...............((((.(((..((.(((((......)))))..))..)))...))))............', -1140),
                        ('...((((.(((..((.(((((......)))))..))..)))...))))...............((((.(((...(.(((((......)))))..)...)))...))))............', -860),
                        ('...((((.(((..((.(((((......)))))..))..)))...))))...............((((.(((.....(((((......)))))......)))...))))............', -970),
                        ('...((((.(((..((.(((((......)))))..))..)))...))))...............((((.((......(((((......))))).......))...))))............', -900),
                        ('...((((.(((..((.(((((......)))))..))..)))...))))...............((((.(.......(((((......)))))........)...))))............', -750),
                        ('...((((.(((..((.(((((......)))))..))..)))...))))...............((((.........(((((......)))))............))))............', -820),
                        ('...((((.(((..((.(((((......)))))..))..)))...))))...............(((((........(((((......)))))........)...))))............', -730),
                        ('...((((.(((..((.(((((......)))))..))..)))...))))...............((((((.......(((((......))))).......))...))))............', -1070),
                        ('...((((.(((..((.(((((......)))))..))..)))...))))...............(((.((.......(((((......))))).......))....)))............', -920),
                        ('...((((.(((..((.(((((......)))))..))..)))...))))...............((((((.......(((((......))))).......)))...)))............', -1180),
                        ('...((((.(((..((.(((((......)))))..))..)))...))))...............(((((((......(((((......)))))......))))...)))............', -1340),
                        ('...((((.(((...(.(((((......)))))..)...)))...))))...............(((((((......(((((......)))))......))))...)))............', -1060), 
                        ('...((((.(((.....(((((......)))))......)))...))))...............(((((((......(((((......)))))......))))...)))............', -1170),
                        ('...((((.((......(((((......))))).......))...))))...............(((((((......(((((......)))))......))))...)))............', -1100),
                        ('...((((.(.......(((((......)))))........)...))))...............(((((((......(((((......)))))......))))...)))............', -950),
                        ('...((((.........(((((......)))))............))))...............(((((((......(((((......)))))......))))...)))............', -1020),
                        ('...(((((........(((((......)))))........)...))))...............(((((((......(((((......)))))......))))...)))............', -930),
                        ('...((((((.......(((((......))))).......))...))))...............(((((((......(((((......)))))......))))...)))............', -1270),
                        ('...(((.((.......(((((......))))).......))....)))...............(((((((......(((((......)))))......))))...)))............', -1120),
                        ('...((((((.......(((((......))))).......)))...)))...............(((((((......(((((......)))))......))))...)))............', -1380),
                        ('...(((((((......(((((......)))))......))))...)))...............(((((((......(((((......)))))......))))...)))............', -1540),
                        ('...(((((((.......((((......)))).......))))...)))...............(((((((......(((((......)))))......))))...)))............', -1420),
                        ('...(((((((.......((((......)))).......))))...)))...............(((((((.......((((......)))).......))))...)))............', -1300)]
        assert barrier == max((en for ss, en in path)) - path[0][1]
        #for (ss, en) in path:
        #    print(f'{ss} {en:>5d}')
        path, barrier = findpath_split(seq, ss1, ss2, vrna_md, th = 5)
        assert barrier == 410
        #for (ss, en) in path:
        #    print(f'{ss} {en:>5d}')
        path, barrier = findpath_split(seq, ss1, ss2, vrna_md, th = 1)
        assert barrier == 410
        #for (ss, en) in path:
        #    print(f'{ss} {en:>5d}')

    def test_findpath_split_02(self):
        seq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCG"
        ss1 = "..((....(((........(((.....)))....)))..(((........)))..))." 
        ss2 = ".((((....))))..(((.(.((...)).)..)))......(((.......)))...."
        #print(f'\n{seq}\n{ss1}\n{ss2}')
        path, barrier = findpath_split(seq, ss1, ss2, RNA.md(), th = 1)
        assert barrier == 200 # may change with a new findpath version?
        #for (ss, en) in path:
        #    print(f'{ss} {en:>5d}')
        #print(barrier)

@unittest.skipIf(SKIP, "skipping tests")
class FloodingTests(unittest.TestCase):
    def test_local_flooding(self):
        seq = "AAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUAC"
        [ss, en] = ['...((........))...........(((((...((((........))))...)))))..', -9.30]
        minh = 1

        # Set model details.
        vrna_md = RNA.md()
        vrna_md.temperature = 25
        fc = RNA.fold_compound(seq, vrna_md)

        macro, fstep = local_flooding(fc, ss, basinh = minh, rates = False)
        for m in macro:
            em = round(fc.eval_structure(m), 2)
            assert em <= en + minh
        for m in fstep:
            em = round(fc.eval_structure(m), 2)
            assert em > en + minh

    def test_path_flooding_barriers(self):
        seq = "AAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUAC"
        path = [('......((.....)).(((((......)))))..((((........))))..........',  -956), # L0013
                ('.......(.....)..(((((......)))))..((((........))))..........',  -587), # S
                ('................(((((......)))))..((((........))))..........',  -897), # I
                ('....(........)..(((((......)))))..((((........))))..........',  -631), # S
                ('...((........)).(((((......)))))..((((........))))..........', -1000), # L0010
                ('...((........))..((((......))))...((((........))))..........',  -890), # I
                ('...((........))...(((......)))....((((........))))..........',  -729), # I
                ('...((........))...((........))....((((........))))..........',  -560), # I
                ('...((........))...(..........)....((((........))))..........',  -288), # S
                ('...((........))...................((((........))))..........',  -681), # I # local minimum
                ('...((........))...............(...((((........))))...)......',  -332), # S
                ('...((........))..............((...((((........))))...)).....',  -517), # I
                ('...((........)).............(((...((((........))))...)))....',  -587), # I
                ('...((........))............((((...((((........))))...))))...',  -699), # I
                ('...((........))...........(((((...((((........))))...)))))..',  -930), # L0017
                ('....(........)............(((((...((((........))))...)))))..',  -561), # S
                ('..........................(((((...((((........))))...)))))..',  -827), # I
                ('.............(....).......(((((...((((........))))...)))))..',  -522), # S
                ('............((....))......(((((...((((........))))...)))))..',  -667), # I
                ('...........(((....))).....(((((...((((........))))...)))))..',  -797), # I
                ('..........((((....))))....(((((...((((........))))...)))))..',  -938), # I
                ('.........(((((....)))))...(((((...((((........))))...)))))..', -1008), # I
                ('.....(...(((((....)))))..)(((((...((((........))))...)))))..',  -820), # S
                ('.....((..(((((....))))).))(((((...((((........))))...)))))..',  -890), # I
                ('.....((...((((....))))..))(((((...((((........))))...)))))..',  -820), # S
                ('.....((.(.((((....))))).))(((((...((((........))))...)))))..',  -969), # I
                ('.....((.(.((((....))))).)).((((...((((........))))...))))...',  -738), # S
                ('....(((.(.((((....))))).)))((((...((((........))))...))))...', -1094)] # L0001

        ssmap = path_flooding(path, minh = 300)

        L0013 = '......((.....)).(((((......)))))..((((........))))..........'
        L0010 = '...((........)).(((((......)))))..((((........))))..........'
        L0017 = '...((........))...........(((((...((((........))))...)))))..'
        L0001 = '....(((.(.((((....))))).)))((((...((((........))))...))))...'
        Lhigh = '...((........))...................((((........))))..........'
        saddles = [1, 8, 10, 17]

        mins = set([L0013, L0010, L0017, L0001, Lhigh])
        seen = set()
        for si in sorted(ssmap):
            lm = ssmap[si]
            if isinstance(lm, list):
                assert si in saddles
                assert len(lm) == 2
                lm1, lm2 = sorted(lm)
                assert path[lm1][0] in mins
                assert path[lm2][0] in mins
            else:
                assert path[lm][0] in mins
                seen.add(path[lm][0])
        assert seen == mins

        ssmap = path_flooding(path, minh = 300, maxlm = -956)
        mins = set([L0013, L0010, L0001])
        seen = set()
        for si in sorted(ssmap):
            lm = ssmap[si]
            if isinstance(lm, list):
                assert si in saddles
                assert len(lm) == 2
                lm1, lm2 = sorted(lm)
            else:
                assert path[lm][0] in mins
                seen.add(path[lm][0])
        assert seen == mins

    def test_path_flooding_random(self):
        seq = "UCUACUAUUCCGGCUUGACAUAAAUAUCGAGUGCUCGACCGCUAUUAUGGUACUUUCCAGCGUUUUGAUUGGUGGAUAAUAUCCCCCAAAAACGCGAGUC"
        path = [('............(((((..........)))))((((..((........)).........((((((...((((.((((...)))).)))))))))))))).', -1850), # lmin
                ('............(((((..........))))).(((..((........)).........((((((...((((.((((...)))).)))))))))))))..', -1710),
                ('............(((((..........)))))..((..((........)).........((((((...((((.((((...)))).))))))))))))...', -1440),
                ('............(((((..........)))))...(..((........)).........((((((...((((.((((...)))).)))))))))))....', -1300), # saddle
                ('............(((((..........)))))......((........)).........((((((...((((.((((...)))).)))))))))).....', -1660), # lmin
                ('............(((((..........))))).......(........)..........((((((...((((.((((...)))).)))))))))).....', -1370),
                ('............(((((..........)))))...........................((((((...((((.((((...)))).)))))))))).....', -1620),
                ('...........((((((..........))))).......)...................((((((...((((.((((...)))).)))))))))).....', -1230),
                ('..........(((((((..........))))).......))..................((((((...((((.((((...)))).)))))))))).....', -1390),
                ('..........((.((((..........))))........))..................((((((...((((.((((...)))).)))))))))).....', -1110), # saddle
                ('..........(((((((..........)))).......)))..................((((((...((((.((((...)))).)))))))))).....', -1520), # lmin
                ('....(.....(((((((..........)))).......)))........).........((((((...((((.((((...)))).)))))))))).....', -1100), # saddle
                ('....((....(((((((..........)))).......))).......)).........((((((...((((.((((...)))).)))))))))).....', -1260),
                ('....(((...(((((((..........)))).......)))......))).........((((((...((((.((((...)))).)))))))))).....', -1380),
                ('....((((..(((((((..........)))).......))).....)))).........((((((...((((.((((...)))).)))))))))).....', -1580),
                ('...(((((..(((((((..........)))).......))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1720), # lmin
                ('...(((((..((((((............))).......))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1440),
                ('...(((((..(((((..............)).......))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1280),
                ('...(((((..((((................).......))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1140), # saddle
                ('...(((((..(((.........................))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1560),
                ('...(((((..(((...(..................)..))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1380),
                ('...(((((..(((..((..................)).))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1480),
                ('...(((((..(((.(((..................)))))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1750),
                ('...(((((..(((.((((................))))))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1870),
                ('...(((((..(((.(((((.............).))))))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1860),
                ('...(((((..(((.((((((...........)).))))))).....)))))........((((((...((((.((((...)))).)))))))))).....', -1950),
                ('...(((((..(((.(((((((.........))).))))))).....)))))........((((((...((((.((((...)))).)))))))))).....', -2150),
                ('..((((((..(((.(((((((.........))).))))))).....)))))).......((((((...((((.((((...)))).)))))))))).....', -2270), # lmin
                ('..((((((..(((.(((((((.........))).))))))).....)))))).....(.((((((...((((.((((...)))).)))))))))))....', -2140),
                ('..((((((..(((.(((((((.........))).))))))).....)))))).....(..(((((...((((.((((...)))).))))))))).)....', -1710)]

        saddles = ['............(((((..........)))))...(..((........)).........((((((...((((.((((...)))).)))))))))))....', 
                   '..........((.((((..........))))........))..................((((((...((((.((((...)))).)))))))))).....',
                   '....(.....(((((((..........)))).......)))........).........((((((...((((.((((...)))).)))))))))).....',
                   '...(((((..((((................).......))).....)))))........((((((...((((.((((...)))).)))))))))).....']

        lmins = ['............(((((..........)))))((((..((........)).........((((((...((((.((((...)))).)))))))))))))).', 
                 '............(((((..........)))))......((........)).........((((((...((((.((((...)))).)))))))))).....',
                 '..........(((((((..........)))).......)))..................((((((...((((.((((...)))).)))))))))).....',
                 '...(((((..(((((((..........)))).......))).....)))))........((((((...((((.((((...)))).)))))))))).....',
                 '..((((((..(((.(((((((.........))).))))))).....)))))).......((((((...((((.((((...)))).)))))))))).....']

        ssmap = path_flooding(path, minh = 300)
        for si in sorted(ssmap):
            lm = ssmap[si]
            if isinstance(lm, list):
                assert len(lm) == 2
                lm1, lm2 = sorted(lm)
                #print(path[si], 'saddle')
                assert path[si][0] in saddles
            elif si == ssmap[si]:
                assert path[si][0] in lmins
                #print(path[si], 'lmin')
            else:
                assert path[si][0] not in lmins
                assert path[si][0] not in saddles
                #print(path[si])
 
        #print() # NOTE: A quick check if edge_flooding works.
        #[s1, e1] = path[0]
        #[s2, e2] = path[-1]
        #fp = init_findpath_max(seq)
        #for (ss1, en1, ssB, enB, ss2, en2) in edge_flooding(fp, s1, s2, e1, e2, minh = 300):
        #    print(ss1, en1, ssB, enB, ss2, en2)

@unittest.skipIf(SKIP, "skipping tests")
class NeighborhoodTests(unittest.TestCase):
    def test_bpd_i_cache(self):
           p = '.((((.((((.((...(((...).))...))))))..))))...'
           i = '.((((.((((........(...)........))))..))))...'
           q = '.((((.((((.((...))(...).((...))))))..))))...'
           assert get_bpd_cache(p, i) == get_bpd_i_cache(p, q)
           assert get_bpd_cache(i, q) == get_bpd_i_cache(q, p)
           assert get_bpd_i_cache(p, i) == get_bpd_i_cache(p, q)
           assert get_bpd_i_cache(q, i) == get_bpd_i_cache(q, p)
           assert get_bpd_i_cache(i, p) == 0
           assert get_bpd_i_cache(i, q) == 0

    def test_guiding_edge_search_01(self):
        sss = ["..........((((.....((((.((.........)).)))).))))...",
               "...((((...)))).....((((.((.........)).))))........",
               ".(((......)))......((((.((.........)).))))........"]
        edges = guiding_edge_search(sss)
        assert len(edges) == 4
        #assert len(edges) == 6

    def test_guiding_edge_search_02(self):
        # NOTE: I did not actually check if this is correct.
        btree = """
              AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGGUGACAACAUU
            1 ..........((((.((((.((.((.......)).))))))..))))...  -6.70    0  13.00
            2 ..........((((.((((.((...((.....)).))))))..))))...  -6.10    1   2.10
            3 ..........((((.....((((.((.........)).)))).))))...  -5.90    1   6.30
            4 ((((.....(((........)))....))))....(((...)))......  -5.70    1   8.50
            5 ...((((...)))).....((((.((.........)).))))........  -5.60    3   4.80
            6 .(((......)))......((((.((.........)).))))........  -5.50    5   4.30
            7 ..........((((..((((.....((.....)).....))))))))...  -5.50    3   5.30
            8 ((((.....(((........)))....))))...................  -5.00    4   3.40
            9 ((((.....((.((.....))))....))))....(((...)))......  -5.00    4   2.80
           10 ((((.....((.((.....)).))...))))....(((...)))......  -4.90    4   3.50
        """
        lmins = parse_barriers(btree, is_file = False, return_tuple = True)
        sti = {x.structure: x.id for x in lmins[1:]}
        edges = guiding_edge_search(sti.keys())
        #print()
        #for (x, y) in sorted(edges, key=lambda x: (sti[x[0]], sti[x[1]])):
        #    if sti[x] < sti[y]:
        #        print(sti[x], sti[y], x, y)
        assert len(edges) == 22
        #assert len(edges) == 24

    def test_guiding_edge_search_03(self):
        seq = "AGACGACAAGGUUGAAUCGCA"
        sss = """(.((......)))........
                 .((((((...))))..))...
                 .(((......)))........
                 .((...((....))..))...
                 .(.(((((....))..)))).
                 .(.((((...))))..)....
                 .(.(((..........)))).
                 .(.((............))).
                 ...(((((....))..)))..
                 ...((((...)))).......
                 ...((((......)..)))..
                 ...(((..........)))..
                 ...((............))..
                 ....................."""
        #sss = """(.((......)))........
        #         ...((((......)..)))..
        #         ...(((..........)))..
        #         ...((............))..
        #         ....................."""
 
        sti = {s: e for e, s in enumerate(sss.split(), 1)}
        #print(len(sti))
        #for ss in sti:
        #    print(sti[ss], ss)
        edges = guiding_edge_search(set(sti.keys()))
        #print()
        #for (x, y) in sorted(edges, key=lambda x: (sti[x[0]], sti[x[1]])):
        #    if sti[x] < sti[y]:
        #        print(sti[x], sti[y], x, y)
        assert len(edges) == 36
        #assert len(edges) == 56

    def test_guide_graph_construction(self):
        btree = """
              AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGGUGACAACAUU
            #1 ..........((((.((((.((.((.......)).))))))..))))...  -6.70    0  13.00
            2 ..........((((.((((.((...((.....)).))))))..))))...  -6.10    1   2.10
            3 ..........((((.....((((.((.........)).)))).))))...  -5.90    1   6.30
            #4 ((((.....(((........)))....))))....(((...)))......  -5.70    1   8.50
            5 ...((((...)))).....((((.((.........)).))))........  -5.60    3   4.80
            6 .(((......)))......((((.((.........)).))))........  -5.50    5   4.30
            #7 ..........((((..((((.....((.....)).....))))))))...  -5.50    3   5.30
            8 ((((.....(((........)))....))))...................  -5.00    4   3.40
            #9 ((((.....((.((.....))))....))))....(((...)))......  -5.00    4   2.80
           #10 ((((.....((.((.....)).))...))))....(((...)))......  -4.90    4   3.50
        """
        lmins = parse_barriers(btree, is_file = False, return_tuple = True)
        seq, md = lmins[0], RNA.md()
        ndata = {x.structure: {'energy': int(round(x.energy*100)), 'identity': x.id} for x in lmins[1:]}

        nodes, edges = get_guide_graph(seq, md, ndata.keys())
        assert all(n not in ndata for n in nodes)
        assert len(edges) == 16

    def test_edge_flooding(self):
        btree = """
              AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGGUGACAACAUU
            1 ..........((((.((((.((.((.......)).))))))..))))...  -6.70    0  13.00
            2 ..........((((.((((.((...((.....)).))))))..))))...  -6.10    1   2.10
            3 ..........((((.....((((.((.........)).)))).))))...  -5.90    1   6.30
            4 ((((.....(((........)))....))))....(((...)))......  -5.70    1   8.50
            5 ...((((...)))).....((((.((.........)).))))........  -5.60    3   4.80
            6 .(((......)))......((((.((.........)).))))........  -5.50    5   4.30
            7 ..........((((..((((.....((.....)).....))))))))...  -5.50    3   5.30
            8 ((((.....(((........)))....))))...................  -5.00    4   3.40
            9 ((((.....((.((.....))))....))))....(((...)))......  -5.00    4   2.80
           10 ((((.....((.((.....)).))...))))....(((...)))......  -4.90    4   3.50
        """
        lmins = parse_barriers(btree, is_file = False, return_tuple = True)
        seq, md = lmins[0], RNA.md()
        sss = [x.structure for x in lmins[1:]]
        sti = {x.structure: x.id for x in lmins[1:]}

        s1 = lmins[1].structure
        s2 = lmins[4].structure
        e1 = int(round(lmins[1].energy*100))
        e2 = int(round(lmins[4].energy*100))

        fp = init_findpath_max(seq)
        #for (ss1, en1, ssB, enB, ss2, en2) in edge_flooding(fp, s1, s2, e1, e2, minh = None):
        #    print(ss1, en1, ssB, enB, ss2, en2)
        assert len(list(edge_flooding((seq, md), s2, s1, e2, e1, minh = None))) == 1

        #for (ss1, en1, ssB, enB, ss2, en2) in edge_flooding(fp, s1, s2, e1, e2, minh = 0):
        #    print(ss1, en1, ssB, enB, ss2, en2)
        assert len(list(edge_flooding(fp, s2, s1, e2, e1, minh = 0))) == 1 # Used to be 6

        #for (ss1, en1, ssB, enB, ss2, en2) in edge_flooding(seq, md, s1, s2, e1, e2, minh = 300):
        #    print(ss1, en1, ssB, enB, ss2, en2)
        assert len(list(edge_flooding(fp, s2, s1, e2, e1, minh = 300))) == 1 # Used to be 3

    def test_guided_neighborhood_flooding_01(self):
        seq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCGUACUUAUACCGCCUGUGCGGACUA"
        sss = """(((((........)))))((((((((((((((((((...(((........))).)))))........)))))))).))))).
                 (((((........)))))(((((((((((((..((....)).....(((((.......)))))....)))))))).))))).
                 .........((((((.((((..((((.........))))..)))).))).))).(((((((...........)))))))...
                 .((((.....(((....)))))))((((((((((((...(((........))).)))))........)))))))........
                 .((((....)))).((((......((((...(((((...(((........))).))))).......)))).....))))...
                 (((((........)))))..((.(((...(((.((....)).)))...))))).(((((((...........)))))))...
                 .((((....)))).......((.(((...(((.((....)).)))...))))).(((((((...........)))))))...
                 (((((........)))))(((((.((((...(((((...(((........))).))))).......))))......))))).
                 (((((........))))).((((...)))).(((((...(((........))).))))).......((((....))))....
                 .((((....)))).((((.((((...))))..))))(((...))).(((((.......)))))...((((....))))....
                 .((((.....(((....)))))))((((...(((((...(((........))).))))).......))))............
                 ........((((..((((.......))))....))))..(((........))).(((((((...........)))))))...
                 .((((...((((..((((.......))))....))))..))))...(((((.......)))))...((((....))))....
                 .((((....)))).((((.......))))..(((((...(((........))).))))).......................
                 ........((((..((((.......))))....))))..(((.((...((((.((........)).))))...)).)))...
                 ......((((((..((((.......))))...................((((.((........)).))))......))))))"""

        fc = RNA.fold_compound(seq, RNA.md())
        ndata = {s: {'identity' : e, 'energy': int(round(fc.eval_structure(s)*100))} for e, s in enumerate(sss.split(), 1)}
        #print(len(ndata))
        #for ss in ndata:
        #    print(ndata[ss], ss)
        gedges = guiding_edge_search(set(ndata.keys()))
        #print()
        #for (x, y) in sorted(gedges, key=lambda x: (ndata[x[0]]['identity'], ndata[x[1]]['identity'])):
        #    if ndata[x]['identity'] < ndata[y]['identity']:
        #        print(ndata[x], ndata[y], x, y)
        edata = dict()
        ndata, edata = neighborhood_flooding((seq, RNA.md()), ndata, gedges, tedges = edata, minh = 300)
        #for (x, y) in sorted(edata, key=lambda x: (ndata[x[0]]['energy'], ndata[x[1]]['energy'])):
        #    print(x, y, edata[(x,y)]['saddle_energy'])
        assert len(gedges) == 44
        assert len(ndata) == 31
        assert len(edata) == 100

    def test_guided_neighborhood_flooding_02(self):
        btree = """
              AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGGUGACAACAUU
            1 ..........((((.((((.((.((.......)).))))))..))))...  -6.70    0  13.00
            2 ..........((((.((((.((...((.....)).))))))..))))...  -6.10    1   2.10
            3 ..........((((.....((((.((.........)).)))).))))...  -5.90    1   6.30
            4 ((((.....(((........)))....))))....(((...)))......  -5.70    1   8.50
            5 ...((((...)))).....((((.((.........)).))))........  -5.60    3   4.80
            6 .(((......)))......((((.((.........)).))))........  -5.50    5   4.30
            7 ..........((((..((((.....((.....)).....))))))))...  -5.50    3   5.30
            8 ((((.....(((........)))....))))...................  -5.00    4   3.40
            9 ((((.....((.((.....))))....))))....(((...)))......  -5.00    4   2.80
           10 ((((.....((.((.....)).))...))))....(((...)))......  -4.90    4   3.50
        """
        lmins = parse_barriers(btree, is_file = False, return_tuple = True)
        seq, md = lmins[0], RNA.md()
        ndata = {x.structure: {'energy': int(round(x.energy*100)), 'identity': x.id} for x in lmins[1:]}
        gnodes, gedges = get_guide_graph(seq, md, ndata.keys())
        assert len(gnodes) == 0
        assert len(gedges) == 36

        fp = init_findpath_max(seq)
        #print()
        #for node in ndata:
        #    print(node, ndata[node])
        nd = {k:v for k,v in ndata.items()}
        ndata, edata = neighborhood_flooding(fp, nd, gedges, minh = 200)
        assert len(ndata) == 11
        assert len(edata) == 40
        #print()
        #for node in ndata:
        #    print(node, ndata[node])

        nd = {k:v for k,v in ndata.items()}
        ndata, edata = neighborhood_flooding(fp, nd, gedges, minh = 100)
        assert len(ndata) == 11
        assert len(edata) == 42

        nd = {k:v for k,v in ndata.items()}
        ndata, edata = neighborhood_flooding(fp, nd, gedges, minh = 500)
        assert len(ndata) == 11
        assert len(edata) == 40

    def test_top_down_coarse_graining_direct_paths(self):
        pass

    def test_top_down_coarse_graining(self):
        seq = "GCCCUUGUCGAGAGGA"
        ss  = "..((((.....))))."

        # Set model details.
        vrna_md = RNA.md()
        fc = RNA.fold_compound(seq, vrna_md)

        macro, fstep = local_flooding(fc, ss, basinh = 9.11, rates = False)
        ndata = {m: {'identity': e,
            'energy': int(round(fc.eval_structure(m)*100))} for e, m in enumerate(sorted(macro, key = lambda x: macro[x][0]))}
        gedges = guiding_edge_search(set(ndata.keys()))
        edata = dict()
        for (x, y) in gedges:
            if get_bpd_cache(x, y) == 1:
                edata[(x,y)] = {'saddle_energy': max(ndata[x]['energy'], ndata[y]['energy'])}
            else:
                edata[(x,y)] = {'saddle_energy': None}

        #assert len(ndata) == 230
        #assert len(edata) == 982
        print()
        print('done flooding', len(ndata), len(edata))

        cgn, cge, cgm = top_down_coarse_graining(ndata, edata, minh = 200)

        print()
        for n in sorted(cgn, key = lambda x: cgn[x]['energy']):
            print(f"{cgn[n]['identity']:>2d} {n} {cgn[n]['energy']:>5d} {len(cgm[n]):>5d}")

        print('Total nodes in mapping', sum(len(cgm[n]) for n in cgn))
        assert sum(len(cgm[n]) for n in cgn) >= len(ndata)-len(cgn)

        for node in sorted(cgn, key = lambda x: cgn[x]['energy']):
            ne = cgn[node]['energy']
            # Calculate barrier heights to all other basins.
            barstr = ''
            for other in sorted(cgn, key = lambda x: cgn[x]['energy']):
                oe = cgn[other]['energy']
                sE = cge[(node, other)]['saddle_energy'] if (node, other) in cge else None
                if sE is not None:
                    barstr += ' {:7.2f}'.format((sE - ne)/100)
                else:
                    barstr += ' {:7.2f}'.format(float('nan'))
            print(barstr)

    def test_minitrafo_randseq(self):
        # A random set of sequences returned by randseq -l 100 | RNAsubopt --stochBT_en=50 | sort -u
        btree = """ 
           CAAAGGCGACUCUCCUUAGACUCUAUAAAUAGUAAAUAGCUCCUAGGGACAAGGCUUACGUCCGCGUUAUUCACAUAAGUCGUCGCUUCAGUUGUGCGAC
        01 ............((((((((..((((.........)))).)).)))))).............................(((((.((.......))))))) -10.40 0 1
        02 .....((((((.(((((((...((.............))...)))))))...(((....)))...((.....))...))))))................. -10.80 0 1
        03 ....((((((..((((((((..((((.........)))).)).))))))...((.(.(((....))).).))......))))))(((......).))... -11.20 0 1
        04 ....(((((((.(((((((...((((.........))))...)))))))..(.((........)).)..........)))))))..(.(.((...)))). -11.60 0 1
        05 .....(((((..(((((((...(((((.....)..))))...)))))))...((((((.((..(.....)..)).))))))))))).............. -13.80 0 1
        06 ...(((((((..((((((....(((...........)))....))))))...((((((.((...........)).)))))))))))))((....)).... -14.20 0 1
        07 ...(.(((.....((((((...((((.........))))...))))))....(((....)))))).)...............((((.........)))). -14.60 0 1
        08 ....((((((...((((((...((((.........))))...))))))....((((((.((...........)).))))))))))))..(....)..... -15.00 0 1
        09 ..((.(((.(..(((((((...((((.........))))...)))))))...)((....)).))).)).............(((((.........))))) -15.10 0 1
        10 ...((....))..((((((...((((.........))))...))))))....((((((.................))))))(((((.........))))) -15.80 0 1
        11 ....((((((..(((((((...((((.........))))...)))))))....((........))((.....)).......))))))............. -16.00 0 1
        12 ....((((((..(((((((...((((.........))))...)))))))...((((((.(.............).))))))))))))...((......)) -16.10 0 1
        13 ....(((((((.(((((((...((.............))...))))))).((.((........)).)).........)))))))((.........))... -16.20 0 1
        14 ..........(.(((((((...((((.........))))...))))))).).((((((.(.............).))))))(((((.........))))) -16.20 0 1
        15 ...((.(.....(((((((...((((.........))))...)))))))...).))......................((((.(((.......))))))) -16.60 0 1
        16 ...(((((((...((((((...((((..(...)..))))...))))))....((((((.((...........)).)))))))))))))............ -17.00 0 1
        17 ....((((((...((((((...((((.........))))...))))))....((((((.((...........)).))))))))))))((.((...)))). -17.00 0 1
        18 ((((((((((..(((((((...((((........).)))...)))))))...((((((.((...........)).)))))))))))))...)))...... -17.20 0 1
        19 ...(((((((..(((((((...((((.........))))...)))))))...((((((.((....(.....))).)))))))))))))((....)).... -17.30 0 1
        20 ...(((((((...(((((((..((((.........)))).).))))))....((((((.((...........)).)))))))))))))............ -17.40 0 1
        21 ((((((((((..(((((((...((((.........))))...)))))))...((((((.((..(......).)).)))))))))))))...)))...... -17.70 0 1
        22 ....((((((..((((((((..((((.........)))).).)))))))...((((((.((...........)).))))))))))))............. -17.90 0 1
        23 .....(((((...((((((...((((.........))))...))))))....((((((.((...........)).))))))))))).............. -17.90 0 1
        24 (..(((((((..(((((((...(((...........)))...)))))))...((((((.((...........)).)))))))))))))..)......... -18.10 0 1
        25 ....(((((((.(((((((...((((.........))))...)))))))..(.((........)).)..........)))))))((.........))... -18.20 0 1
        26 (..(((((((..(((((((...((((.........))))...)))))))...((((((((....))).........))))))))))))..)......... -18.30 0 1
        27 ...(((((((..((((((((..((((.........)))).).)))))))...((((((.((...........)).)))))))))))))............ -18.50 0 1
        28 ((((((((((..(((((((...((((.........))))...)))))))...((((((.(.............).)))))))))))))...)))...... -18.70 0 1
        29 ...((((.....(((((((...((((.........))))...)))))))....))))........................(((((.........))))) -18.70 0 1
        30 ....((((((..(((((((...((((.........))))...)))))))...((((((.((...........)).))))))))))))...(((....))) -18.70 0 1
        31 ...(((((((..(((((((...((((.........))))...)))))))...((((((.((...........)).)))))))))))))..((...))... -18.90 0 1
        32 ....((((((..(((((((...((((.........))))...)))))))...((((((.((...........)).))))))))))))((........)). -18.90 0 1
        33 .....((((((.(((((((...((((.........))))...))))))).((.((........)).)).........))))))(((.........))).. -18.90 0 1
        34 ....((((((...((((((...((((.........))))...))))))....((((((.((...........)).))))))))))))............. -19.20 0 1
        35 ..((.(((....(((((((...((((.........))))...)))))))...(((....)))))).)).............(((((.........))))) -19.30 0 1
        36 ...(((((((...((((((...((((.........))))...))))))....((((((.((...........)).)))))))))))))............ -19.80 0 1
        37 ....((((((..(((((((...((((.........))))...)))))))...((((((.((...........)).))))))))))))............. -20.30 0 1
        38 ...(((((((..(((((((...((((.........))))...)))))))...((((((.((...........)).)))))))))))))............ -20.90 0 1
        39 ....(((((((.(((((((...((((.........))))...)))))))....((........))............)))))))................ -17.60 0 1
        """
        lmins = parse_barriers(btree, is_file = False, return_tuple = True)
        seq, md = lmins[0], RNA.md()
        myminh = 300

        #print()
        ndata = {x.structure: {'energy': int(round(x.energy*100)), 'identity': x.id} for x in lmins[1:]}

        #print(f'Finding guide neighborhood for {len(ndata)=}.')
        gnodes, gedges = get_guide_graph(seq, md, ndata.keys())
        #print(f' - Found {len(gedges)} guide edges and {len(gnodes)} new guide nodes.')
        assert len(gnodes) == 25
        assert len(gedges) == 240
        for nid, (ss, en) in enumerate(gnodes, 40):
            ndata[ss] = {'energy': en, 'identity': 40}

        fp = init_findpath_max(seq)
        ndata, edata = neighborhood_flooding(fp, ndata, gedges, minh = myminh)

        # Those results are not constant ... why?
        #assert len(ndata) == 86
        #assert len(edata) == 340

        print(len(ndata), len([(x,y) for (x, y) in edata if edata[(x,y)]['saddle_energy'] is not None]))
        cgn, cge, cgm = top_down_coarse_graining(ndata, edata, minh = myminh)

        print()
        for n in sorted(cgn, key = lambda x: cgn[x]['energy']):
            print(f" {n} {cgn[n]['energy']:>5d} {len(cgm[n]):>5d}")

        print('Total nodes in mapping', sum(len(cgm[n]) for n in cgn))
        for node in sorted(cgn, key = lambda x: cgn[x]['energy']):
            ne = cgn[node]['energy']
            # Calculate barrier heights to all other basins.
            barstr = ''
            for other in sorted(cgn, key = lambda x: cgn[x]['energy']):
                oe = cgn[other]['energy']
                sE = cge[(node, other)]['saddle_energy'] if (node, other) in cge else None
                if sE is not None:
                    barstr += ' {:7.2f}'.format((sE - ne)/100)
                else:
                    barstr += ' {:7.2f}'.format(float('nan'))
            print(barstr)


        ##print()
        ##for n in sorted(ndata, key = lambda x: ndata[x]['energy']):
        ##    print(n, ndata[n]['energy'])
        ##print()
        ##for n in sorted(cg_ndata, key = lambda x: cg_ndata[x]['energy']):
        ##    print(n, cg_ndata[n]['energy'], len(cg_ndata[n]['hiddennodes']))

        #hiddennodes = set()
        #for n in cg_ndata:
        #    if cg_ndata[n]['hiddennodes']:
        #        hiddennodes |= cg_ndata[n]['hiddennodes']
        #assert len(cg_ndata) + len(hiddennodes) == len(ndata)




if __name__ == '__main__':
    unittest.main()
