#!/usr/bin/env python
import os
import math
import shutil
import tempfile
import unittest
import networkx as nx

import RNA
from ribolands.syswraps import sys_treekin_051
from ribolands.parser import parse_barriers
import ribolands.trafo as trafo
import ribolands.utils as ru

skip = False

def write_log(TL, nlist):
    tlen = len(TL.transcript)
    for e, (ni, data) in enumerate(nlist, 1):
        print("{:4d} {:4d} {} {:6.2f} {:6.4f} ID = {:d}".format(tlen, e, ni[:tlen], 
                data['energy'], data['occupancy'], data['identity']))

@unittest.skipIf(skip, "slow tests are disabled by default")
class Test_TrafoLand(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix='DrTest_')
        print("\nWriting temporary files to: {}".format(self.tmpdir))

    def tearDown(self):
        print("Removing temporary file directory: {}".format(self.tmpdir))
        #shutil.rmtree(self.tmpdir)

    def _init_TL(self, seq, sss, add_edges = True):
        fullseq = seq
        vrna_md = RNA.md()

        TL = trafo.TrafoLandscape(fullseq, vrna_md)
        TL._transcript_length = len(seq)

        for ss in sss:
            if not TL.has_node(ss):
                TL.addnode(ss, structure = ss)
                TL.nodes[ss]['occupancy'] = 1.0 / len(sss),

        return TL

    def test_edge_attributes(self):
        """testing:
        TrafoLandscape.has_edge()
        TrafoLandscape.get_rate()
        TrafoLandscape.get_saddle()
        """
        fullseq = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCA"
        TL = trafo.TrafoLandscape(fullseq, RNA.md())
        s1 = ".....(((((......)).)))((((((((....)))))....)))..."
        s2 = "..................((((.(((((((....)))))....))))))"

        TL.add_weighted_edges_from([(s1, s2, 0)])
        TL[s1][s2]['saddleE'] = float('inf')

        self.assertTrue(TL.has_edge(s1, s2))
        self.assertEqual(TL.get_saddle(s1, s2), float('inf'))
        self.assertEqual(TL.get_rate(s1, s2), 0)

        TL[s1][s2]['saddleE'] = 9999
        TL[s1][s2]['weight'] = 0.9999
        self.assertEqual(TL.get_saddle(s1, s2), 9999)
        self.assertEqual(TL.get_rate(s1, s2), 0.9999)

    def test_TL_fraying_helices(self):
        seq = "CCCCGGAAGGAAUGGGUAGUGACAGCAGCUUAGGUCGCUGCAUCAUCCCC"
        sss = """
        .............(((...((...(((((.......)))))..))..)))
        ....(((.....(((((((((((..........)))))))).))))))..
        ........((.((((((((((((..........)))))))).))))))..
        ........((.((((((((((((..(......))))))))).)))).)).
        ........((.(((.((((((((..........))))))))..))).)).
        ........(((.(((((((((((..........)))))))).))))))..
        ........(((.(((((((((((((...))...)))))))).))))))..
        .....(..(((....((((((((..(......)))))))))....))).)
        ........(((.((.((((((((((...))...))))))))..)))))..
        ..((....))...(((..((((..(((((.......))))).))))))).
        ..((....))...((((.(((.((((.((....)).))))))).))))..
        ..((....))...(((..((((((((.((....)).))))..))))))).
        ........((...((((.(((.((((.((....)).))))))).))))))
        ((......))...((((.(((.((((.((....)).))))))).).))).
        """.split()

        TL = self._init_TL(seq, sss)
        TL.minh = 3
        snodes = TL.sorted_nodes(attribute = 'energy', rev = True)
        enth = TL.nodes[snodes[0]]['energy']

        nn = TL.connect_nodes_n2()
        TL.coarse_grain()
        print(len(TL.active_nodes))

        mss, mfe = TL.fc.backtrack(len(seq))
        if TL.has_node(mss):
            assert TL.nodes[mss]['active']

        nn = TL.connect_nodes_n2(nodes = TL.active_nodes)
        print(nn)

        TL.coarse_grain()
        print(len(TL.active_nodes), len(TL.nodes))

        nbrs = TL.find_fraying_neighbors(TL.active_nodes, mfree = 6)
        for par, new in sorted(nbrs.items(), key = lambda x: x[0]):
            TL.connect_nodes_n2(nodes = new)
        TL.coarse_grain()
        print(len(TL.active_nodes), len(TL.nodes))

    def test_TL_connect_mfe(self):
        rbar = """# A random barriers example:
             CCCCGGAAGGAAUGGGUAGUGACAGCAGCUUAGGUCGCUGCAUCAUCCCC
           1 ........((.((((((((((((..........)))))))).)))).)). -15.10    0  15.00
           2 ..((....))...((((.(((.((((.((....)).))))))).)))).. -15.00    1  13.50
           3 ..((....))...(((..((((..(((((.......))))).))))))). -14.90    2  11.90
           4 ........(((.(((((((((((..........)))))))).)))))).. -14.90    1   3.30
           5 ..((....))...(((..((((..(((((.......))))).)))).))) -14.40    3   5.60
           6 ..((....))...(((..(((.((((.((....)).)))))))...))). -14.30    2   6.20
           7 ..((....)).((((((((((((..........)))))))).)))).... -13.70    1   3.30
           8 .((.....))...((((.(((.((((.((....)).))))))).)))).. -13.60    2   3.30
           9 .((((.......))))..((((..(((((.......))))).)))).... -13.50    3   7.20
          10 .((.....))...(((..((((..(((((.......))))).))))))). -13.50    3   3.30
          11 ..((....))...(((..((((((((.((....)).))))..))))))). -13.50    6   4.20
          12 .((.....))...(((..((((..(((((.......))))).)))).))) -13.00    5   3.30
          13 ..((....))...(((..((((((((.((....)).))))..)))).))) -13.00    2   5.60
          14 .((.....))...(((..(((.((((.((....)).)))))))...))). -12.90    6   3.30
        """

        lmins = parse_barriers(rbar, is_file = False, return_tuple = True)

        seq = lmins[0]
        mss = lmins[1].structure
        mfe = float(lmins[1].energy)

        sss = [] # Everything except for MFE structure
        for l in lmins[2:]:
            sss.append(l.structure)

        TL = self._init_TL(seq, sss)
        TL.minh = 3.30
        for n in TL.nodes():
            TL.nodes[n]['active'] = True

        assert mss not in sss
        assert not TL.has_node(mss)

        nn = TL.connect_to_active(nodes = [mss])
        for node in nn:
            assert TL.nodes[node]['active'] is None

        assert TL.nodes[mss]['active'] is None
        assert nx.is_strongly_connected(TL)

    def test_minitrafo(self, verbose = False):
        fullseq = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCA"
        TL = trafo.TrafoLandscape(fullseq, RNA.md())
        fname = self.tmpdir + '/' + 'minitrafo'

        self.assertEqual(list(TL.nodes), [])
        self.assertEqual(TL.transcript, '')
        self.assertEqual(TL._transcript_length, 0)
        TL.expand()
        self.assertEqual(len(TL), 1)
        self.assertEqual(list(TL.nodes), ['.' * len(fullseq)])
        self.assertEqual(TL.nodeID, 1)

        [bfile, rfile, p0, nlist] = TL.get_simulation_files_tkn(fname)

        self.assertEqual(TL._transcript_length, 1)
        self.assertEqual(TL.transcript, 'C')
        stepsize = 2

        for i in range(2, len(fullseq), stepsize):
            seq  = fullseq[0:i]
            tlen = len(TL.transcript)
            TL.expand(extend = i - tlen)
            TL.coarse_grain()
            self.assertEqual(i, len(TL.transcript))
            [bfile, rfile, p0, nlist] = TL.get_simulation_files_tkn(fname)
            if len(nlist) == 1:
                TL.total_time += 0.2
            else:
                bfile = None # sometimes bfile causes a segfault, ...
                tfile = sys_treekin(fname, rfile, 
                        treekin = 'treekin', 
                        bofile = bfile,
                        binrates = True, 
                        p0 = p0, t0 = 0, ti = 1.5, t8 = 0.2, 
                        exponent = False, 
                        useplusI = False, 
                        force = True)

                time_inc, iterations = TL.update_occupancies_tkn(tfile, nlist)
                TL.total_time += time_inc

            if verbose:
                write_log(TL, nlist)

            dn, sr = TL.prune(0.01)

    def dont_test_expand_and_coarse_grain(self, verbose = True):
        seq = "AUAUAGCUUGUUUACUUUGGAUGAACUGGGGAGAAAAUCCUGGUAAAACU"
        sss = [
            "..........((((((..((((...((....))...)))).))))))...",
            ".......(((((((...)))))))((((((........))))))......",
            ".......(((((((...)))))))((((((.......))).)))......",
            ".(((......(((.((((((.....))))))..)))......)))....."
        ]

        TL = self._init_TL(seq, sss)

        ess = [
            ["..........((((((..((((...((....))...)))).))))))...", 0.25, True],
            [".......(((((((...)))))))((((((........))))))......", 0.25, True],
            [".......(((((((...)))))))((((((.......))).)))......", 0.25, True],
            [".(((......(((.((((((.....))))))..)))......))).....", 0.25, True]]
        for (ss, occ, active) in ess:
            self.assertTrue(TL.has_node(ss))
            self.assertEqual(TL.nodes[ss]['occupancy'], occ)
            self.assertEqual(TL.nodes[ss]['active'], active)

        TL.expand(extend=0, exp_mode='fullconnect')
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
            self.assertTrue(TL.has_node(ss))
            self.assertEqual(TL.nodes[ss]['occupancy'], occ)
            self.assertEqual(TL.nodes[ss]['active'], active)

        if verbose:
            TL.get_simulation_files_tkn(self.tmpdir+'/ecp1')

        TL.coarse_grain(dG_min=4.3)
        ess = [
            ["..........((((((..((((...((....))...)))).))))))...", 0.25, True],
            ["........................((((((........))))))......", 0.125, True],
            [".......(((((((...)))))))((((((........))))))......", 0.375, True],
            [".........((((((....).))))).(((.......)))..........", 0.00, True],
            ["........................((((((.......))).)))......", 0.00, False],
            [".......(((((((...)))))))((((((.......))).)))......", 0.00, False],
            [".....(((..(((.((((((.....))))))..))).....)))......", 0.25, True],
            [".(((......(((.((((((.....))))))..)))......))).....", 0.00, False]]
        for (ss, occ, active) in ess:
            print(ss, occ, active)
            print(TL.nodes[ss])
            self.assertTrue(TL.has_node(ss))
            self.assertEqual(TL.nodes[ss]['occupancy'], occ)
            self.assertEqual(TL.nodes[ss]['active'], active)
        
        if verbose:
            TL.get_simulation_files_tkn(self.tmpdir+'/ecp2')

        TL.prune(0.01)
        ess = [
            ["..........((((((..((((...((....))...)))).))))))...", 0.25, True],
            ["........................((((((........))))))......", 0.125, True],
            [".......(((((((...)))))))((((((........))))))......", 0.375, True],
            [".........((((((....).))))).(((.......)))..........", 0.00, False],
            ["........................((((((.......))).)))......", 0.00, False],
            [".......(((((((...)))))))((((((.......))).)))......", 0.00, False],
            [".....(((..(((.((((((.....))))))..))).....)))......", 0.25, True],
            [".(((......(((.((((((.....))))))..)))......))).....", 0.00, False]]
        for (ss, occ, active) in ess:
            print(ss, occ, active)
            self.assertTrue(TL.has_node(ss))
            self.assertEqual(TL.nodes[ss]['occupancy'], occ)
            self.assertEqual(TL.nodes[ss]['active'], active)

        if verbose:
            TL.get_simulation_files_tkn(self.tmpdir+'/ecp3')

    def dont_test_coarse_graining_dG(self, verbose = False):
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
        bar = parse_barriers(bfile, return_tuple = True)

        fullseq = "UGAAUGUGCCGCUAGACGACAUCCCGCCGGAUGGCGGGGC"
        vrna_md = RNA.md()

        CG = trafo.TrafoLandscape(fullseq, vrna_md)
        CG._transcript_length = len(fullseq)

        min_dG = 3.0

        bmap = dict()
        inactive = set() 
        for lmin in bar[1:]:
            id = lmin.id
            ss = lmin.structure
            en = lmin.energy
            fa = lmin.father
            ba = lmin.barrier

            bmap[ss]=id
            if fa:
                ss2 = bar[fa][1]
                en2 = round(bar[fa].energy, 2)
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
            CG.nodes[ss]['active'] = True
            CG.nodes[ss2]['active'] = True
            CG.nodes[ss]['energy'] = float(en)
            CG.nodes[ss2]['energy'] = float(en2)
            CG.nodes[ss]['occupancy'] = 0
            CG.nodes[ss2]['occupancy'] = 0

        for lmin in bar[1:]:
            id = lmin.id
            ss = lmin.structure
            en = lmin.energy
            fa = lmin.father
            ba = lmin.barrier

            if verbose:
                print('a', id, ss, en, fa, ba, CG.nodes[ss]['active'])
            nbrs = filter(lambda x: CG.nodes[x]['active'], sorted(CG.successors(ss), 
                              key=lambda x: (CG.nodes[x]['energy'], x), reverse=False))
            self.assertTrue(CG.has_node(ss))
            self.assertEqual(CG.nodes[ss]['active'], True)

        mn = CG.coarse_grain(dG_min=min_dG)

        if verbose:
            for x in sorted(mn, key=lambda x:int(bmap[x])):
                sm = '-> {}'.format([bmap[y] for y in mn[x]])
                print('b', "{} {} {:s}".format(x, bmap[x], sm))

        for lmin in bar[1:]:
            id = lmin.id
            ss = lmin.structure
            en = lmin.energy
            fa = lmin.father
            ba = lmin.barrier
            if verbose:
                print('c', id, ss, en, fa, ba, CG.nodes[ss]['active'])
            self.assertTrue(CG.has_node(ss))
            if ss in inactive:
                self.assertFalse(CG.nodes[ss]['active'])
            else:
                self.assertEqual(CG.nodes[ss]['active'], True)

    def dont_test_coarse_graining_by_rates(self, verbose = False):
        """
        A coarse graining test based on this example...
        All structures are connected based on rates in the corresponding 
        rates file. Because coarse-graining requires a saddle energy, it is
        estimated using dG = -RT * log(rate), but that is not consistent 
        with the barrier heights found in the bar file...
        Particularly, 3 gets merged to 2, at dG = 2.9, because 4 gets merged
        into 3, ...

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
        bar = parse_barriers(bfile, return_tuple = True)
        rts = ru.parse_ratefile(rfile)

        fullseq = "UGAAUGUGCCGCUAGACGACAUCCCGCCGGAUGGCGGGGC"
        vrna_md = RNA.md()

        CG = trafo.TrafoLandscape(fullseq, vrna_md)
        CG._transcript_length = len(fullseq)

        if verbose:
            for x in bar:
                print(x)

        min_dG = 4

        inactive = set() 
        for m1, line in enumerate(rts, 1):
            s1 = bar[m1].structure
            en1 = round(bar[m1].energy, 2)
            for m2, rate in enumerate(line, 1):
                if m1 == m2 :
                    continue
                if rate :
                    s2 = bar[m2].structure
                    en2 = round(bar[m2].energy, 2)

                    CG.add_weighted_edges_from([(s1, s2, rate)])

                    # this is an approximation, the barriers rate is not computed 
                    # from a single transition state....
                    dG = -CG._RT * math.log(rate)
                    sE = round(en1+dG, 2)
                    CG[s1][s2]['saddle'] = sE

                    if en1 > en2 and dG < min_dG:
                        inactive.add(s1)
                        if verbose:
                            print('{} -> {} | {:.2f} {:.5f} {} {}'.format(
                                m1+1, m2+1, dG, rate, bar[m1].barrier, sE))

            CG.nodes[s1]['active'] = True
            CG.nodes[s1]['energy'] = en1
            CG.nodes[s1]['occupancy'] = 0

        for lmin in bar[1:]:
            id = lmin.id
            ss = lmin.structure
            en = lmin.energy
            fa = lmin.father
            ba = lmin.barrier
            #print('b', id, ss, en, fa, ba, CG.nodes[ss]['active'])
            nbrs = filter(lambda x: CG.nodes[x]['active'], sorted(CG.successors(ss), 
                              key=lambda x: (CG.nodes[x]['energy'], x), reverse=False))
            #for (x,y) in zip(nbrs, map(lambda x: CG[ss][x]['saddle'], nbrs)):
            #    print('   ', x, y)
            self.assertTrue(CG.has_node(ss))
            self.assertTrue(CG.nodes[ss]['active'])

        CG.coarse_grain(dG_min=min_dG)

        for lmin in bar[1:]:
            id = lmin.id
            ss = lmin.structure
            en = lmin.energy
            fa = lmin.father
            ba = lmin.barrier
 
            nbrs = filter(lambda x: CG.nodes[x]['active'], sorted(CG.successors(ss), 
                              key=lambda x: (CG.nodes[x]['energy'], x), reverse=False))
            if verbose:
                print(id, ss, en, fa, ba, CG.nodes[ss]['active'], end='')
                print(' ', list(map(lambda x: CG[ss][x]['saddle'], nbrs)))
            self.assertTrue(CG.has_node(ss))
            if ss in inactive:
                self.assertFalse(CG.nodes[ss]['active'])
            elif int(id) == 3:
                # This is interesting, staring from min_dG=2.9 you can see that 
                # the assertTrue would break.. because based on transition rates,
                # there exists a path 3->4->2 which is energetically favorable 
                # over 3->2, so after 4 is merged into 3, 3 can be merged into 2.
                self.assertFalse(CG.nodes[ss]['active'])
            else:
                self.assertTrue(CG.nodes[ss]['active'])

@unittest.skipIf(skip, "slow tests are disabled by default")
class Test_HelperFunctions(unittest.TestCase):
    def test_fold_exterior_loop(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        ss = ".....(((((......)).)))((((((((....))))....)))).."

        ext_moves = dict()
        md = RNA.md()

        nbr, ext = trafo.fold_exterior_loop(se, ss, md, ext_moves)

        self.assertEqual(nbr, '.....(((((......)).)))((((((((....))))....))))..')
        self.assertEqual(ext, 'CUCGUCGNNNCGGGNNNCCGC')
        self.assertEqual(ext_moves[ext], '.....((...))((...))..')

    def test_fold_exterior_loop_2(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        ss = ".....(((((......)).)))((((((((....))))....)))).."

        md = RNA.md()
        nbr, ext = trafo.fold_exterior_loop(se, ss, md)

        self.assertEqual(nbr, '.....(((((......)).)))((((((((....))))....))))..')
        self.assertEqual(ext, 'CUCGUCGNNNCGGGNNNCCGC')
        self.assertEqual(trafo.EFOLD_CACHE[ext], '.....((...))((...))..')

    def test_open_fraying_helices(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        ss = ".....(((((......)).)))((((((((....))))....)))).."

        out = trafo.open_fraying_helices(se, ss, free=6)

        res = [
            "........((......))....((((((((....))))....))))..",
            ".....(((((......)).)))....((((....))))..........",
            "........((......))........((((....)))).........."]

        self.assertEqual(sorted(out), sorted(res))

        out = trafo.open_fraying_helices(se, ss, free=8)

        res = [
            "......................((((((((....))))....))))..",
            ".....(((((......)).)))....((((....))))..........",
            "..........................((((....)))).........."]

        self.assertEqual(sorted(out), sorted(res))

    def test_open_fraying_helices_multi(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCG"
        ss = "..((.(((((......)).)))((((((((....))))....)))).))"

        out = trafo.open_fraying_helices(se, ss, free=6)
        res = [".....(((((......)).)))((((((((....))))....))))..."]
        self.assertEqual(sorted(out), sorted(res))

        out = trafo.open_fraying_helices(se, ss, free=7)
        res = [
            "........((......))....((((((((....))))....))))...",
            ".....(((((......)).)))....((((....))))...........",
            "........((......))........((((....))))..........."]
        self.assertEqual(sorted(out), sorted(res))

if __name__ == '__main__':
    unittest.main()
