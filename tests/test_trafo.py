#!/usr/bin/env python
import math
import shutil
import tempfile
import unittest

import RNA
from ribolands.parser import parse_barriers
from ribolands.utils import parse_ratefile

import ribolands.trafo as trafo
from ribolands.trafo import (TrafoLandscape,
                             fold_exterior_loop, 
                             open_fraying_helices)

SKIP = False

def write_log(TL, nlist):
    tlen = len(TL.transcript)
    for e, (ni, data) in enumerate(nlist, 1):
        print("{:4d} {:4d} {} {:6.2f} {:6.4f} ID = {:d}".format(tlen, e, ni[:tlen], 
                data['energy'], data['occupancy'], data['identity']))

@unittest.skipIf(SKIP, "slow tests are disabled by default")
class Test_TrafoLand(unittest.TestCase):
    def _init_TL(self, seq, sss):
        fullseq = seq
        vrna_md = RNA.md()

        TL = TrafoLandscape(fullseq, vrna_md)
        TL._transcript_length = len(seq)

        for ss in sss:
            if not TL.has_node(ss):
                TL.addnode(ss, structure = ss, occupancy = 1/len(sss))
        return TL

    def test_edge_attributes(self):
        """testing:
        TrafoLandscape.has_edge()
        TrafoLandscape.get_rate()
        TrafoLandscape.get_saddle()
        """
        fullseq = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCA"
        TL = TrafoLandscape(fullseq, RNA.md())
        s1 = ".....(((((......)).)))((((((((....)))))....)))..."
        s2 = "..................((((.(((((((....)))))....))))))"
        
        TL.addnode(s1, structure = s1)
        TL.addnode(s2, structure = s2)
        TL.addedge(s1, s2, weight = 0)
        TL.edges[(s1, s2)]['weight'] == 0

        self.assertTrue(TL.has_edge(s1, s2))
        self.assertEqual(TL.get_rate(s1, s2), 0)

        TL.edges[(s1, s2)]['weight'] = 0.9999
        self.assertEqual(TL.get_rate(s1, s2), 0.9999)

    def test_TL_expansion_random(self):
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
        TL.minh = 300
        mss, mfe = TL.fc.backtrack(len(seq))

        #print()
        #for sn in TL.sorted_nodes(attribute = 'energy'):
        #    print(sn, TL.nodes[sn]['energy'])

        nn = TL.expand()
        #print()
        #for sn in TL.sorted_nodes(attribute = 'energy'):
        #    print(sn, TL.nodes[sn]['energy'])
        assert len(nn) + len(sss) == len(TL.nodes)
        assert mss in TL.nodes
        
    def test_TL_expansion_btree(self):
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
        TL.minh = 330
        TL.maxh = 1180 # lower leads to disconnected components
        assert not TL.has_node(mss)

        #print()
        #for sn in TL.sorted_nodes(attribute = 'energy'):
        #    print(sn, TL.nodes[sn]['energy'])
        assert len(list(TL.sccs())) == 13

        nn = TL.expand()
        assert len(nn) + len(sss) == len(TL.nodes)
        assert mss in TL.nodes

        #print()
        #for sn in TL.sorted_nodes(attribute = 'energy'):
        #    print(sn, TL.nodes[sn]['energy'])
        
        #for ed in sorted(TL.edges):
        #    print(ed, TL.edges[ed]['saddle_energy'])
        # TODO: did not check those results
        assert len(list(TL.sccs())) == 1
        assert len(list(TL.sccs(minh =  330))) == 14
        assert len(list(TL.sccs(minh =  500))) ==  7
        assert len(list(TL.sccs(minh =  800))) ==  3
        assert len(list(TL.sccs(minh = 1000))) ==  3
        assert len(list(TL.sccs(minh = 1150))) ==  3
        assert len(list(TL.sccs(minh = 1550))) ==  1

    def test_expand_and_coarse_grain(self, verbose = True):
        seq = "AUAUAGCUUGUUUACUUUGGAUGAACUGGGGAGAAAAUCCUGGUAAAACU"
        sss =["..........((((((..((((...((....))...)))).))))))...",
              ".......(((((((...)))))))((((((........))))))......",
              ".......(((((((...)))))))((((((.......))).)))......",
              ".(((......(((.((((((.....))))))..)))......)))....."]

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

        TL.expand(extend=0)

        ess = [
            ['..........((((((..((((...((....))...)))).))))))...', 0.25, True], # -860 0.25
            ['........................((((((........))))))......', 0.00, True], # -590 0
            ['.......(((((((...)))))))((((((........))))))......', 0.25, True], # -570 0.25
            ['..................((((...((....))...))))..........', 0.00, True], # -430 0
            ['.........((((((....).))))).(((.......)))..........', 0.00, True], # -230 0
            ['........................((((((.......))).)))......', 0.00, True], # -220 0
            ['.......(((((((...)))))))((((((.......))).)))......', 0.25, True], # -200 0.25
            ['...........................(((........))).........', 0.00, True], # -170 0
            ['.....(((..(((.((((((.....))))))..))).....)))......', 0.00, True], # -140 0
            ['..............((((((.....))))))...................', 0.00, True], # -140 0
            ['.(((......(((.((((((.....))))))..)))......))).....', 0.25, True]] # -70  0.25

        for (ss, occ, active) in ess:
            self.assertTrue(TL.has_node(ss))
            self.assertEqual(TL.nodes[ss]['occupancy'], occ)
            self.assertEqual(TL.nodes[ss]['active'], active)

        cn, ce = TL.get_coarse_network(minh = 430)
        #print()
        #for n in sorted(TL.nodes, key = lambda x: TL.nodes[x]['energy']):
        #    print(n, TL.nodes[n]['energy'], TL.nodes[n]['occupancy'], TL.nodes[n]['active'])

        ess = [
            ['..........((((((..((((...((....))...)))).))))))...', 0.25, True], # -860 0.25
            ['........................((((((........))))))......', 0.00, True], # -590 0
            ['.......(((((((...)))))))((((((........))))))......', 0.25, True], # -570 0.25
            ['..................((((...((....))...))))..........', 0.00, False], # -430 0
            ['.........((((((....).))))).(((.......)))..........', 0.00, True], # -230 0
            ['........................((((((.......))).)))......', 0.00, True], # -220 0
            ['.......(((((((...)))))))((((((.......))).)))......', 0.25, True], # -200 0.25
            ['...........................(((........))).........', 0.00, False], # -170 0
            ['.....(((..(((.((((((.....))))))..))).....)))......', 0.25, True], # -140 0
            ['..............((((((.....))))))...................', 0.00, False], # -140 0
            ['.(((......(((.((((((.....))))))..)))......))).....', 0.00, False]] # -70  0.25

        for (ss, occ, active) in ess:
            self.assertTrue(TL.has_node(ss))
            self.assertEqual(TL.nodes[ss]['occupancy'], occ)
            self.assertEqual(TL.nodes[ss]['active'], active)
        
        TL.prune(0.01)

        ess = [
            ['..........((((((..((((...((....))...)))).))))))...', 0.25, True], # -860 0.25
            ['........................((((((........))))))......', 0.00, False], # -590 0
            ['.......(((((((...)))))))((((((........))))))......', 0.25, True], # -570 0.25
            ['..................((((...((....))...))))..........', 0.00, False], # -430 0
            ['.........((((((....).))))).(((.......)))..........', 0.00, False], # -230 0
            ['........................((((((.......))).)))......', 0.00, False], # -220 0
            ['.......(((((((...)))))))((((((.......))).)))......', 0.25, True], # -200 0.25
            ['...........................(((........))).........', 0.00, False], # -170 0
            ['.....(((..(((.((((((.....))))))..))).....)))......', 0.25, True], # -140 0
            ['..............((((((.....))))))...................', 0.00, False], # -140 0
            ['.(((......(((.((((((.....))))))..)))......))).....', 0.00, False]] # -70  0.25

        for (ss, occ, active) in ess:
            self.assertTrue(TL.has_node(ss))
            self.assertEqual(TL.nodes[ss]['occupancy'], occ)
            self.assertEqual(TL.nodes[ss]['active'], active)

@unittest.skipIf(SKIP, "slow tests are disabled by default")
class Test_HelperFunctions(unittest.TestCase):
    def test_fold_exterior_loop(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        ss = ".....(((((......)).)))((((((((....))))....)))).."

        ext_moves = dict()
        md = RNA.md()

        nbr, ext = fold_exterior_loop(se, ss, md, ext_moves)

        self.assertEqual(nbr, '.....(((((......)).)))((((((((....))))....))))..')
        self.assertEqual(ext, 'CUCGUCGNNNCGGGNNNCCGC')
        self.assertEqual(ext_moves[ext], '.....((...))((...))..')

    def test_fold_exterior_loop_2(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        ss = ".....(((((......)).)))((((((((....))))....)))).."

        md = RNA.md()
        nbr, ext = fold_exterior_loop(se, ss, md)

        self.assertEqual(nbr, '.....(((((......)).)))((((((((....))))....))))..')
        self.assertEqual(ext, 'CUCGUCGNNNCGGGNNNCCGC')
        self.assertEqual(trafo.EFOLD_CACHE[ext], '.....((...))((...))..')

    def test_open_fraying_helices(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        ss = ".....(((((......)).)))((((((((....))))....)))).."

        out = open_fraying_helices(se, ss, free=6)

        res = [
            "........((......))....((((((((....))))....))))..",
            ".....(((((......)).)))....((((....))))..........",
            "........((......))........((((....)))).........."]

        self.assertEqual(sorted(out), sorted(res))

        out = open_fraying_helices(se, ss, free=8)

        res = [
            "......................((((((((....))))....))))..",
            ".....(((((......)).)))....((((....))))..........",
            "..........................((((....)))).........."]

        self.assertEqual(sorted(out), sorted(res))

    def test_open_fraying_helices_multi(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCG"
        ss = "..((.(((((......)).)))((((((((....))))....)))).))"

        out = open_fraying_helices(se, ss, free=6)
        res = [".....(((((......)).)))((((((((....))))....))))..."]
        self.assertEqual(sorted(out), sorted(res))

        out = open_fraying_helices(se, ss, free=7)
        res = [
            "........((......))....((((((((....))))....))))...",
            ".....(((((......)).)))....((((....))))...........",
            "........((......))........((((....))))..........."]
        self.assertEqual(sorted(out), sorted(res))

if __name__ == '__main__':
    unittest.main()
