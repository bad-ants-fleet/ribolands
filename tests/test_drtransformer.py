#!/usr/bin/env python

import sys
import unittest
import math

import RNA
import ribolands.trafo as trafo
from ribolands.syswraps import sys_treekin
import ribolands.utils as ru


class Test_TrafoLandscape(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def dont_test_minitrafo(self):
        # remove the dont_ for testing, but beware it writes files ...
        fullseq = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCA"
        vrna_md = RNA.md()

        CG = trafo.TrafoLandscape(fullseq, vrna_md)

        self.assertEqual(list(CG.nodes), [])

        CG._transcript_length = 0
        CG.expand()
        self.assertEqual(len(CG), 1)
        self.assertEqual(list(CG.nodes), ['.' * len(fullseq)])
        self.assertEqual(CG._nodeid, 1)

        [bfile, rfile, p0, nlist] = CG.get_simulation_files_tkn('rudi')

        self.assertEqual(CG._transcript_length, 1)

        stepsize = 2

        for i in range(2, len(fullseq), stepsize):
            seq = fullseq[0:i]
            #print i, CG._transcript_length, i-CG._transcript_length
            CG.expand(extend=i-CG._transcript_length)
            self.assertEqual(i, CG._transcript_length)
            # if i == 10:
            #  CG.logfile = sys.stdout
            [bfile, rfile, p0, nlist] = CG.get_simulation_files_tkn('rudi')
            if len(nlist) == 1:
                CG._total_time += 0.2
            else:
                # sometimes bfile causes a segfault, so let's leave it out.
                bfile = None
                tfile, _ = sys_treekin('rudi', seq, bfile, rfile, binrates=True,
                                       treekin='treekin', p0=p0, t0=0, ti=1.5, t8=0.2,
                                       exponent=False, useplusI=False, force=True, verb=False)

                time_inc, iterations = CG.update_occupancies_tkn(tfile, nlist)
                CG._total_time += time_inc

            dn, sr, rj = CG.prune()

    def test_interface(self):
        pass

    def test_transition_edge(self):
        fullseq = "GGAACCGUCUCCCUCUGCCAAAAGGUAGAGGGAGAUGGAGCAUCUCUCUCUACGAAGCAGAGAGAGACGAAGG"
        vrna_md = RNA.md()

        CG = trafo.TrafoLandscape(fullseq, vrna_md)

        self.assertEqual(list(CG.nodes()), [])

        CG.expand()
        self.assertEqual(len(CG), 1)
        self.assertEqual(list(CG.nodes()), ['.' * len(fullseq)])
        self.assertEqual(CG._nodeid, 1)

    def test_expand_and_coarse_grain(self):

        seq = "AUAUAGCUUGUUUACUUUGGAUGAACUGGGGAGAAAAUCCUGGUAAAACU"
        sss = [
            "..........((((((..((((...((....))...)))).))))))...",
            ".......(((((((...)))))))((((((........))))))......",
            ".......(((((((...)))))))((((((.......))).)))......",
            ".(((......(((.((((((.....))))))..)))......)))....."
        ]

        CG = self.initialize_CG(seq, sss)

        ess = [
            ["..........((((((..((((...((....))...)))).))))))...", 0.25, True],
            [".......(((((((...)))))))((((((........))))))......", 0.25, True],
            [".......(((((((...)))))))((((((.......))).)))......", 0.25, True],
            [".(((......(((.((((((.....))))))..)))......))).....", 0.25, True]]
        for (ss, occ, active) in ess:
            self.assertTrue(CG.has_node(ss))
            self.assertEqual(CG.node[ss]['occupancy'], occ)
            self.assertEqual(CG.node[ss]['active'], active)

        CG.expand(extend=0)
        ess = [
            ["..........((((((..((((...((....))...)))).))))))...", 0.25, True],
            ["........................((((((........))))))......", 0.00, True],
            [".......(((((((...)))))))((((((........))))))......", 0.25, True],
            [".........((((((....).))))).(((.......)))..........", 0.00, True],
            ["........................((((((.......))).)))......", 0.00, True],
            [".......(((((((...)))))))((((((.......))).)))......", 0.25, True],
            [".....(((..(((.((((((.....))))))..))).....)))......", 0.00, True],
            [".(((......(((.((((((.....))))))..)))......))).....", 0.25, True]]
        for (ss, occ, active) in ess:
            self.assertTrue(CG.has_node(ss))
            self.assertEqual(CG.node[ss]['occupancy'], occ)
            self.assertEqual(CG.node[ss]['active'], active)

        # CG.get_simulation_files_tkn('expand')

        CG.coarse_grain(dG_min=4.3)
        ess = [
            ["..........((((((..((((...((....))...)))).))))))...", 0.25, True],
            ["........................((((((........))))))......", 0.00, True],
            [".......(((((((...)))))))((((((........))))))......", 0.50, True],
            [".........((((((....).))))).(((.......)))..........", 0.00, True],
            ["........................((((((.......))).)))......",
             0.00, False],
            [".......(((((((...)))))))((((((.......))).)))......", 0.00, False],
            [".....(((..(((.((((((.....))))))..))).....)))......", 0.25, True],
            [".(((......(((.((((((.....))))))..)))......))).....", 0.00, False]]
        for (ss, occ, active) in ess:
            self.assertTrue(CG.has_node(ss))
            self.assertEqual(CG.node[ss]['occupancy'], occ)
            self.assertEqual(CG.node[ss]['active'], active)

        CG.prune()
        #CG.logfile = sys.stdout
        # CG.get_simulation_files_tkn('expand_again')
        ess = [
            ["..........((((((..((((...((....))...)))).))))))...", 0.25, True],
            #["........................((((((........))))))......", 0.00, False],
            ["........................((((((........))))))......", 0.00, True],
            [".......(((((((...)))))))((((((........))))))......", 0.50, True],
            [".........((((((....).))))).(((.......)))..........", 0.00, False],
            ["........................((((((.......))).)))......", 0.00, False],
            [".......(((((((...)))))))((((((.......))).)))......", 0.00, False],
            [".....(((..(((.((((((.....))))))..))).....)))......", 0.25, True],
            [".(((......(((.((((((.....))))))..)))......))).....", 0.00, False]]
        for (ss, occ, active) in ess:
            self.assertTrue(CG.has_node(ss))
            self.assertEqual(CG.node[ss]['occupancy'], occ)
            self.assertEqual(CG.node[ss]['active'], active)

    def initialize_CG(self, seq, sss):
        fullseq = seq
        vrna_md = RNA.md()

        CG = trafo.TrafoLandscape(fullseq, vrna_md)
        CG._transcript_length = len(seq)

        for e, s1 in enumerate(sss, 1):
            if not CG.has_node(s1):
                en = round(CG._fold_compound.eval_structure(s1), 2)
                CG.add_node(s1, energy=en, occupancy=1.0 / len(sss),
                            identity=CG._nodeid, active=True, last_seen=0)
                CG._nodeid += 1
            for s2 in sss[e:]:
                if not CG.has_node(s2):
                    en = round(CG._fold_compound.eval_structure(s2), 2)
                    CG.add_node(s2, energy=en, occupancy=1.0 / len(sss),
                                identity=CG._nodeid, active=True, last_seen=0)
                    CG._nodeid += 1
                assert CG.add_transition_edge(s1, s2)

        return CG

    def test_coarse_graining_dG(self):
        """
        A coarse graining test based on this example...
        All structures are connected based on the barrier-tree, 
        using Arrhenius rates.

             UGAAUGUGCCGCUAGACGACAUCCCGCCGGAUGGCGGGGC
           1 .....((.((((((..((.(.....).))..)))))).)) -14.90    0  20.00
           2 .....(((.((.....)).)))((((((....)))))).. -13.90    1  13.60
           3 ...((((...........))))((((((....)))))).. -13.60    2   3.40
           4 ...((((..((.....))))))((((((....)))))).. -13.20    3   2.10
           5 .....................(((((((....))))))). -13.20    3   2.80
           6 .....((........))....(((((((....))))))). -12.90    2   3.60
           7 ...(((((.(....).)).)))((((((....)))))).. -12.80    5   1.90
           8 .....((.((((.......(((((....))))))))).)) -12.80    1   6.50
           9 .....((.(((((.......((((....))))))))).)) -11.80    8   2.60
          10 .....((.((((((.......(((....))))))))).)) -11.70    9   2.20
          11 .....((.((((((..((......)).....)))))).)) -11.00    1   2.50
          12 .....((...)).........(((((((....))))))). -11.00    2   1.20
          13 ..........((..........((((((....)))))))) -10.40    2   2.60
          14 .......((((.....))....((((((....)))))))) -10.20    2   2.30
          15 ((..((......))..))...(((((((....))))))).  -9.90    2   1.10
          16 .......(((.............(((((....))))))))  -9.80    2   4.00
          17 .......((.............((((((....))))))))  -9.80    2   1.90
          18 ....(.((....)).).....(((((((....))))))).  -9.60    2   1.20
          19 .((.(((...........)))))(((((....)))))...  -9.50    2   2.70
          20 .......(((....((.....))(((((....))))))))  -9.30   16   1.30
        """

        bfile = 'tests/files/ex1.bar'
        bar = ru.parse_barfile(bfile)

        fullseq = "UGAAUGUGCCGCUAGACGACAUCCCGCCGGAUGGCGGGGC"
        vrna_md = RNA.md()

        CG = trafo.TrafoLandscape(fullseq, vrna_md)
        CG._transcript_length = len(fullseq)

        min_dG = 3.0

        bmap = dict()
        inactive = set() 
        for (id, ss, en, fa, ba) in bar:
            bmap[ss]=id
            if int(fa):
                ss2 = bar[int(fa)-1][1]
                en2 = round(float(bar[int(fa)-1][2]),2)
            else :
                continue
            
            # ensure saddle is not lower than s1, s2
            saddleE = round(float(en) + float(ba),2)

            # Energy barrier
            dG_1s = saddleE - float(en)
            dG_2s = saddleE - float(en2)

            if dG_1s < min_dG+0.001:
                inactive.add(ss)

            # Metropolis Rule
            k_12 = CG._k0 * math.exp(-dG_1s / CG._RT)
            k_21 = CG._k0 * math.exp(-dG_2s / CG._RT)
            CG.add_weighted_edges_from([(ss, ss2, k_12)])
            CG.add_weighted_edges_from([(ss2, ss, k_21)])
            CG[ss][ss2]['saddle'] = saddleE
            CG[ss2][ss]['saddle'] = saddleE
            CG.node[ss]['active'] = True
            CG.node[ss2]['active'] = True
            CG.node[ss]['energy'] = float(en)
            CG.node[ss2]['energy'] = float(en2)
            CG.node[ss]['occupancy'] = 0
            CG.node[ss2]['occupancy'] = 0

        for (id, ss, en, fa, ba) in bar:
            #print 'a', id, ss, en, fa, ba, CG.node[ss]['active']
            #nbrs = filter(lambda x: CG.node[x]['active'], sorted(CG.successors(ss), 
            #                  key=lambda x: (CG.node[x]['energy'], x), reverse=False))
            #for (x,y) in zip(nbrs, map(lambda x: CG[ss][x]['saddle'], nbrs)):
            #    print '   ', x, y
            self.assertTrue(CG.has_node(ss))
            self.assertEqual(CG.node[ss]['active'], True)

        mn = CG.coarse_grain(dG_min=min_dG)

        #for x in sorted(mn, key=lambda x:int(bmap[x])):
        #    print bmap[x], bmap[mn[x]]

        for (id, ss, en, fa, ba) in bar:
            #print id, ss, en, fa, ba, CG.node[ss]['active']
            self.assertTrue(CG.has_node(ss))
            if ss in inactive:
                self.assertFalse(CG.node[ss]['active'])
            else:
                self.assertEqual(CG.node[ss]['active'], True)



    def test_coarse_graining_by_rates(self):
        """
        A coarse graining test based on this example...
        All structures are connected based on rates in the corresponding 
        rates file. Because coarse-graining requires a saddle energy, it is
        estimated using dG = -RT * log(rate), but that is not consistent 
        with the barrier heights found in the bar file...
        Particularly, 3 gets merged to 2, at dG = 2.9, because 4 gets merged
        into 3, and then also 

             UGAAUGUGCCGCUAGACGACAUCCCGCCGGAUGGCGGGGC
           1 .....((.((((((..((.(.....).))..)))))).)) -14.90    0  20.00
           2 .....(((.((.....)).)))((((((....)))))).. -13.90    1  13.60
           3 ...((((...........))))((((((....)))))).. -13.60    2   3.40
           4 ...((((..((.....))))))((((((....)))))).. -13.20    3   2.10
           5 .....................(((((((....))))))). -13.20    3   2.80
           6 .....((........))....(((((((....))))))). -12.90    2   3.60
           7 ...(((((.(....).)).)))((((((....)))))).. -12.80    5   1.90
           8 .....((.((((.......(((((....))))))))).)) -12.80    1   6.50
           9 .....((.(((((.......((((....))))))))).)) -11.80    8   2.60
          10 .....((.((((((.......(((....))))))))).)) -11.70    9   2.20
          11 .....((.((((((..((......)).....)))))).)) -11.00    1   2.50
          12 .....((...)).........(((((((....))))))). -11.00    2   1.20
          13 ..........((..........((((((....)))))))) -10.40    2   2.60
          14 .......((((.....))....((((((....)))))))) -10.20    2   2.30
          15 ((..((......))..))...(((((((....))))))).  -9.90    2   1.10
          16 .......(((.............(((((....))))))))  -9.80    2   4.00
          17 .......((.............((((((....))))))))  -9.80    2   1.90
          18 ....(.((....)).).....(((((((....))))))).  -9.60    2   1.20
          19 .((.(((...........)))))(((((....)))))...  -9.50    2   2.70
          20 .......(((....((.....))(((((....))))))))  -9.30   16   1.30
        """

        bfile = 'tests/files/ex1.bar'
        rfile = 'tests/files/ex1.rts'
        bar = ru.parse_barfile(bfile)
        rts = ru.parse_ratefile(rfile)

        fullseq = "UGAAUGUGCCGCUAGACGACAUCCCGCCGGAUGGCGGGGC"
        vrna_md = RNA.md()

        CG = trafo.TrafoLandscape(fullseq, vrna_md)
        CG._transcript_length = len(fullseq)

        #for (id, ss, en, fa, ba) in bar:
        #    print 'aa', id, ss, en, fa, ba

        min_dG = 4

        inactive = set() 
        for m1, line in enumerate(rts):
            s1 = bar[m1][1]
            en1 = round(float(bar[m1][2]), 2)
            for m2, rate in enumerate(line):
                if m1 == m2 :
                    continue
                if rate :
                    s2 = bar[m2][1]
                    en2 = round(float(bar[m2][2]), 2)

                    CG.add_weighted_edges_from([(s1, s2, rate)])

                    # this is an approximation, the barriers rate is not computed 
                    # from a single transition state....
                    dG = -CG._RT * math.log(rate)
                    sE = round(en1+dG, 2)
                    CG[s1][s2]['saddle'] = sE

                    if en1 > en2 and dG < min_dG:
                        inactive.add(s1)
                        #print m1+1, m2+1, s1, s2, dG, rate, bar[m1][4], sE

            CG.node[s1]['active'] = True
            CG.node[s1]['energy'] = en1
            CG.node[s1]['occupancy'] = 0

        for (id, ss, en, fa, ba) in bar:
            #print 'a', id, ss, en, fa, ba, CG.node[ss]['active']
            #nbrs = filter(lambda x: CG.node[x]['active'], sorted(CG.successors(ss), 
            #                  key=lambda x: (CG.node[x]['energy'], x), reverse=False))
            #for (x,y) in zip(nbrs, map(lambda x: CG[ss][x]['saddle'], nbrs)):
            #    print '   ', x, y
            self.assertTrue(CG.has_node(ss))
            self.assertEqual(CG.node[ss]['active'], True)

        CG.coarse_grain(dG_min=min_dG)

        for (id, ss, en, fa, ba) in bar:
            #print id, ss, en, fa, ba, CG.node[ss]['active']
            nbrs = filter(lambda x: CG.node[x]['active'], sorted(CG.successors(ss), 
                              key=lambda x: (CG.node[x]['energy'], x), reverse=False))
            #print nbrs
            self.assertTrue(CG.has_node(ss))
            if ss in inactive:
                self.assertFalse(CG.node[ss]['active'])
            else:
                # This is interesting, staring from min_dG=2.9 you can see that 
                # the assertTrue would break.. because based on transition rates,
                # there exists a path 3->4->2 which is energetically favorable 
                # over 3->2, so after 4 is merged into 3, 3 can be merged into 2.
                #self.assertEqual(CG.node[ss]['active'], True)
                #print map(lambda x: CG[ss][x]['saddle'], nbrs)
                pass

        return


#@unittest.skipIf(True, "slow tests are disabled by default")
class test_stuff(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_fold_exterior_loop(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        ss = ".....(((((......)).)))((((((((....))))....)))).."

        ext_moves = dict()
        md = RNA.md()

        nbr, ext = trafo.fold_exterior_loop(md, se, ss, ext_moves)

        self.assertEqual(
            nbr, '.....(((((......)).)))((((((((....))))....))))..')
        self.assertEqual(ext, 'CUCGUCGNNNCGGGNNNCCGC')
        self.assertEqual(ext_moves[ext][1], '.....((...))((...))..')
        self.assertEqual(ext_moves[ext][0], set())

    def test_open_breathing_helices(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        ss = ".....(((((......)).)))((((((((....))))....)))).."

        out = trafo.open_breathing_helices(se, ss, free=6)

        res = [
            "........((......))....((((((((....))))....))))..",
            ".....(((((......)).)))....((((....))))..........",
            "........((......))........((((....)))).........."]

        self.assertEqual(sorted(out), sorted(res))

        out = trafo.open_breathing_helices(se, ss, free=8)

        res = [
            "......................((((((((....))))....))))..",
            ".....(((((......)).)))....((((....))))..........",
            "..........................((((....)))).........."]

        self.assertEqual(sorted(out), sorted(res))

    def test_open_breathing_helices_multi(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCG"
        ss = "..((.(((((......)).)))((((((((....))))....)))).))"

        out = trafo.open_breathing_helices(se, ss, free=6)
        res = [".....(((((......)).)))((((((((....))))....))))..."]
        self.assertEqual(sorted(out), sorted(res))

        out = trafo.open_breathing_helices(se, ss, free=7)
        res = [
            "........((......))....((((((((....))))....))))...",
            ".....(((((......)).)))....((((....))))...........",
            "........((......))........((((....))))..........."]
        self.assertEqual(sorted(out), sorted(res))

if __name__ == '__main__':
    unittest.main()
