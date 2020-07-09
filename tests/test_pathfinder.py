#!/usr/bin/env python3

import RNA
import unittest

from ribolands.pathfinder import (PathfinderError,
                                  path_flooding, 
                                  local_flooding, 
                                  get_fpath_cache,
                                  get_fpath_flooding_cache,
                                  clear_fpath_cache, 
                                  apply_bp_change,
                                  get_bp_change)


class Test_base_pair_changes(unittest.TestCase):
    def tearDown(self):
        clear_fpath_cache()

    def test_apply_bp_change_bug1(self):
        seq = 'AAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUAC' 
        s1 =  '...((((.......(((((((......)))).......)))...))))............'
        s2 = None
        subseq =    'GCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUU'
        subs1  =    '(......((((((((......)))).......))))..)'
        subs2  =    '((((.......((((......)))).......)))...)'

        # Index error
        with self.assertRaises(IndexError):
            apply_bp_change(seq, s1, s2, subseq, subs1, subs2)

    def test_apply_bp_change_bug2(self):
        seq = 'UAGCAGGAUCGAGUAGACCCGAGGAUUGUAUAGUUAGUCUCCACGGCAAUGACGAUGAGUGUCGAGAAUU'
        s1  = '.....((.((.....)).))..((((((......))))))...((.......))................'
        s2  = None
        subseq = 'AGGNNNCUCCACGGCAAUGACGAUGAG'
        subs1  = '.((xxx))...((.......)).....'
        subs2  = '..(xxx)((..((.......))..)).'
        apply_bp_change(seq, s1, s2, subseq, subs1, subs2)

    def test_apply_bp_change_multiple_results(self):
        # TODO: The routine breaks if there are multiple 
        # changes that have the same bp motif.
        # This is not a feature!
        seq = 'AAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUACAAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUAC'
        s1  = '...((((.(((..((.(((((......)))))..))..)))...))))...............((((.(((..((.(((((......)))))..))..)))...))))............'
        s2  = '....((((((.......((((......)))).......))))...)).................((((((.......((((......)))).......))))...)).............',
        subseq = 'AGCNNNGUA'
        subs1  = '.((xxx)).' 
        subs2  = '..(xxx)..'

        with self.assertRaises(NotImplementedError):
            apply_bp_change(seq, s1, s2, subseq, subs1, subs2)

class Test_Pathflooding(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        clear_fpath_cache()

    def test_get_bp_change(self):
        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        s1 = ".....(((((......)).)))((((((((....))))....)))).."
        s2 = "..................(.....................)......."

        features = get_bp_change(se, s1, s2)
        self.assertEqual(len(features), 1)
        assert features[0][0] == "UCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCG"
        assert features[0][1] == ".(((((......)).)))((((((((....))))....))))."
        assert features[0][2] == "..............(.....................)......"

        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        s1 = ".....(((((......)).)))((((((((....))))....)))).."
        s2 = "........((......))(.......((((....))))..)......."

        features = get_bp_change(se, s1, s2)
        self.assertEqual(len(features), 1)
        assert features[0][0] == "UCGCCNNNGUGCGGGCGCNNNGUUAUCGCCG"
        assert features[0][1] == ".((((xxx).)))(((((xxx)....))))."
        assert features[0][2] == "....(xxx)(.......(xxx)..)......"

        se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
        s1 = ".....(..((......))...)...(((((....))))....)....."
        s2 = "........((......))(.......((((....))))..)......."

        features = get_bp_change(se, s1, s2)
        self.assertEqual(len(features), 1)
        assert features[0][0] == "UCGCCNNNGUGCGGGCGCNNNGUUAUCG"
        assert features[0][1] == ".(..(xxx)...)...((xxx)....)."
        assert features[0][2] == "....(xxx)(.......(xxx)..)..."

    def test_local_flooding(self):
        seq = "AAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUAC"
        [ss, en] = ['...((........))...........(((((...((((........))))...)))))..', -9.30]
        minh = 1

        # Set model details.
        vrna_md = RNA.md()
        vrna_md.noLP = 0
        vrna_md.temperature = 25
        vrna_md.dangles = 2
        vrna_md.logML = 0
        vrna_md.special_hp = 1
        vrna_md.noGU = 0
        
        fc = RNA.fold_compound(seq, vrna_md)

        macro, fstep = local_flooding(fc, ss, basinh = minh, rates = False)

        for m in macro:
            em = round(fc.eval_structure(m), 2)
            assert em <= en + minh
        
        for m in fstep:
            em = round(fc.eval_structure(m), 2)
            assert em > en + minh

    def test_path_flooding(self):
        """testing:"""
        seq = "AAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUAC"
        path = [('......((.....)).(((((......)))))..((((........))))..........',  -9.56), # L0013
                ('.......(.....)..(((((......)))))..((((........))))..........',  -5.87), # S
                ('................(((((......)))))..((((........))))..........',  -8.97), # I
                ('....(........)..(((((......)))))..((((........))))..........',  -6.31), # S
                ('...((........)).(((((......)))))..((((........))))..........', -10.00), # L0010
                ('...((........))..((((......))))...((((........))))..........',  -8.90), # I
                ('...((........))...(((......)))....((((........))))..........',  -7.29), # I
                ('...((........))...((........))....((((........))))..........',  -5.60), # I
                ('...((........))...(..........)....((((........))))..........',  -2.88), # S
                ('...((........))...................((((........))))..........',  -6.81), # I # local minimum
                ('...((........))...............(...((((........))))...)......',  -3.32), # S
                ('...((........))..............((...((((........))))...)).....',  -5.17), # I
                ('...((........)).............(((...((((........))))...)))....',  -5.87), # I
                ('...((........))............((((...((((........))))...))))...',  -6.99), # I
                ('...((........))...........(((((...((((........))))...)))))..',  -9.30), # L0017
                ('....(........)............(((((...((((........))))...)))))..',  -5.61), # S
                ('..........................(((((...((((........))))...)))))..',  -8.27), # I
                ('.............(....).......(((((...((((........))))...)))))..',  -5.22), # S
                ('............((....))......(((((...((((........))))...)))))..',  -6.67), # I
                ('...........(((....))).....(((((...((((........))))...)))))..',  -7.97), # I
                ('..........((((....))))....(((((...((((........))))...)))))..',  -9.38), # I
                ('.........(((((....)))))...(((((...((((........))))...)))))..', -10.08), # I
                ('.....(...(((((....)))))..)(((((...((((........))))...)))))..',  -8.20), # S
                ('.....((..(((((....))))).))(((((...((((........))))...)))))..',  -8.90), # I
                ('.....((...((((....))))..))(((((...((((........))))...)))))..',  -8.20), # S
                ('.....((.(.((((....))))).))(((((...((((........))))...)))))..',  -9.69), # I
                ('.....((.(.((((....))))).)).((((...((((........))))...))))...',  -7.38), # S
                ('....(((.(.((((....))))).)))((((...((((........))))...))))...', -10.94)] # L0001

        lmins, ssmap = path_flooding(path, minh = 3)

        L0013 = '......((.....)).(((((......)))))..((((........))))..........'
        L0010 = '...((........)).(((((......)))))..((((........))))..........'
        L0017 = '...((........))...........(((((...((((........))))...)))))..'
        L0001 = '....(((.(.((((....))))).)))((((...((((........))))...))))...'
        Lhigh = '...((........))...................((((........))))..........'

        assert len(lmins) == 5
        assert L0001 in lmins
        assert L0010 in lmins
        assert L0013 in lmins
        assert L0017 in lmins
        assert Lhigh in lmins

        #for k,v in lmins.items():
        #    print(k, v)

        for k,v in sorted(ssmap.items()):
            #print(k, v)
            if k in [1, 8, 10, 17]:
                assert isinstance(v, list)
            else:
                assert isinstance(v, int)
 

    def test_findpath_flooding_minh(self):
        """testing:
           GCCGCCUUAAGCCUACUUAGAUGGAAGUGACGNNNCACGAUUUU 
           ...((.....)).(((((......)))))..(xxx)........ 
           .(((.(.((((....))))).)))((((...(xxx)...)))). 

         0 ...((.....)).(((((......)))))..(xxx)........   2.30
         1 ....(.....)..(((((......)))))..(xxx)........   5.70
         2 .............(((((......)))))..(xxx)........   2.10
         3 ..............((((......))))...(xxx)........   3.30
         4 ...............(((......)))....(xxx)........   4.70
         5 ................((......)).....(xxx)........   7.40
         6 .................(......)......(xxx)........   8.30
         7 ...............................(xxx)........   3.90
         8 ..........(....)...............(xxx)........   7.40
         9 .........((....))..............(xxx)........   6.10
        10 ........(((....))).............(xxx)........   5.10
        11 .......((((....))))............(xxx)........   3.60
        12 ..(....((((....))))...)........(xxx)........   6.00
        13 .((....((((....))))...)).......(xxx)........   2.90
        14 .(((...((((....))))..))).......(xxx)........   2.10
        15 .(((...((((....))))..)))...(...(xxx)...)....   5.20
        16 .(((...((((....))))..)))..((...(xxx)...))...   3.70
        17 .(((.(.((((....))))).)))..((...(xxx)...))...   2.60
        18 .(((.(.((((....))))).))).(((...(xxx)...)))..   2.00
        19 .(((.(.((((....))))).)))((((...(xxx)...)))).   1.10
        """

        seq = "AAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUAC"
        s1 = '......((.....)).(((((......)))))..((((........))))..........'
        s2 = '....(((.(.((((....))))).)))((((...((((........))))...))))...'


        #  0 ...((.....)).(((((......)))))..(xxx)........   2.30 M
        #  1 ....(.....)..(((((......)))))..(xxx)........   5.70 S (0 -> 3.4)
        #  2 .............(((((......)))))..(xxx)........   2.10 M
        #  3 ..............((((......))))...(xxx)........   3.30 
        #  4 ...............(((......)))....(xxx)........   4.70
        #  5 ................((......)).....(xxx)........   7.40
        #  6 .................(......)......(xxx)........   8.30 S (7 -> 4.4, 2 -> 6.20)
        #  7 ...............................(xxx)........   3.90 M
        #  8 ..........(....)...............(xxx)........   7.40 S (7 -> 3.5)
        #  9 .........((....))..............(xxx)........   6.10
        # 10 ........(((....))).............(xxx)........   5.10
        # 11 .......((((....))))............(xxx)........   3.60 M
        # 12 ..(....((((....))))...)........(xxx)........   6.00 S (11 -> 2.4)
        # 13 .((....((((....))))...)).......(xxx)........   2.90 
        # 14 .(((...((((....))))..))).......(xxx)........   2.10 M
        # 15 .(((...((((....))))..)))...(...(xxx)...)....   5.20 S (14 -> 3.1)
        # 16 .(((...((((....))))..)))..((...(xxx)...))...   3.70
        # 17 .(((.(.((((....))))).)))..((...(xxx)...))...   2.60 OP
        # 18 .(((.(.((((....))))).))).(((...(xxx)...)))..   2.00
        # 19 .(((.(.((((....))))).)))((((...(xxx)...)))).   1.10 M

        def count_steps(seq, s1, s2, md, minh, count):
            lmp = get_fpath_flooding_cache(seq, s1, s2, md, minh)
            if lmp[0] is None:
                if isinstance(lmp[1], set):
                    for (lsq, pm1, pm2) in lmp[1]: 
                        count = count_steps(lsq, pm1, pm2, md, minh, count)
                else:
                    for (pm1, pm2) in lmp[1]: 
                        count = count_steps(seq, pm1, pm2, md, minh, count)
                return count
            else:
                return count + 1

        # All basins separated by a barrier less than minh paramter should be merged
        steps = count_steps(seq, s1, s2, RNA.md(), minh = 2.4, count = 0)
        #self._print_multi_path(seq, s1, s2, seq, s1, s2)
        assert steps == 6 # There is an additional orthogonal step 16->17
        clear_fpath_cache()

        steps = count_steps(seq, s1, s2, RNA.md(), minh = 2.41, count = 0)
        assert steps == 5
        clear_fpath_cache()

        steps = count_steps(seq, s1, s2, RNA.md(), minh = 3.1, count = 0)
        assert steps == 5
        clear_fpath_cache()

        steps = count_steps(seq, s1, s2, RNA.md(), minh = 3.11, count = 0)
        #self._print_multi_path(seq, s1, s2, seq, s1, s2)
        assert steps == 3
        clear_fpath_cache()

        steps = count_steps(seq, s1, s2, RNA.md(), minh = 3.41, count = 0)
        assert steps == 3
        clear_fpath_cache()

        steps = count_steps(seq, s1, s2, RNA.md(), minh = 3.5, count = 0)
        assert steps == 3
        clear_fpath_cache()

        steps = count_steps(seq, s1, s2, RNA.md(), minh = 3.51, count = 0)
        assert steps == 2
        clear_fpath_cache()

        steps = count_steps(seq, s1, s2, RNA.md(), minh = 6.2, count = 0)
        assert steps == 2
        clear_fpath_cache()

        steps = count_steps(seq, s1, s2, RNA.md(), minh = 6.21, count = 0)
        assert steps == 1
        clear_fpath_cache()

    def test_findpath_flooding_nr1(self):
        """testing:
           GCCGCCUUAAGCCUACUUAGAUGGAAGUGACGNNNCACGAUUUU 
           ...((.....)).(((((......)))))..(xxx)........ 
           .(((.(.((((....))))).)))((((...(xxx)...)))). 

         0 ...((.....)).(((((......)))))..(xxx)........   2.30
         1 ....(.....)..(((((......)))))..(xxx)........   5.70
         2 .............(((((......)))))..(xxx)........   2.10
         3 ..............((((......))))...(xxx)........   3.30
         4 ...............(((......)))....(xxx)........   4.70
         5 ................((......)).....(xxx)........   7.40
         6 .................(......)......(xxx)........   8.30
         7 ...............................(xxx)........   3.90
         8 ..........(....)...............(xxx)........   7.40
         9 .........((....))..............(xxx)........   6.10
        10 ........(((....))).............(xxx)........   5.10
        11 .......((((....))))............(xxx)........   3.60
        12 ..(....((((....))))...)........(xxx)........   6.00
        13 .((....((((....))))...)).......(xxx)........   2.90
        14 .(((...((((....))))..))).......(xxx)........   2.10
        15 .(((...((((....))))..)))...(...(xxx)...)....   5.20
        16 .(((...((((....))))..)))..((...(xxx)...))...   3.70
        17 .(((.(.((((....))))).)))..((...(xxx)...))...   2.60
        18 .(((.(.((((....))))).))).(((...(xxx)...)))..   2.00
        19 .(((.(.((((....))))).)))((((...(xxx)...)))).   1.10
        """

        seq = "AAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUAC"
        s1 = '......((.....)).(((((......)))))..((((........))))..........'
        s2 = '....(((.(.((((....))))).)))((((...((((........))))...))))...'

        FLOOD = 3
        md = RNA.md()

        x = get_fpath_flooding_cache(seq, s1, s2, md, minh = FLOOD)
        print(x)

        for (tsq, ts1, ts2) in x[1]:
            fw = get_fpath_flooding_cache(tsq, ts1, ts2, md, minh = FLOOD)
            rv = get_fpath_flooding_cache(tsq, ts2, ts1, md, minh = FLOOD)
            assert fw[0] is None
            assert rv[0] is None
            assert isinstance(fw[1], list)
            assert isinstance(rv[1], list)
            assert fw[1][0][0] == rv[1][-1][1]
            assert fw[1][0][1] == rv[1][-1][0]
            assert fw[1][1][0] == rv[1][-2][1]
            assert fw[1][1][1] == rv[1][-2][0]
            assert fw[1][2][0] == rv[1][-3][1]
            assert fw[1][2][1] == rv[1][-3][0]
            assert fw[1][3][0] == rv[1][-4][1]
            assert fw[1][3][1] == rv[1][-4][0]

        y = get_fpath_flooding_cache(seq, s2, s1, md, minh = FLOOD)
        assert x[0] is None
        assert y[0] is None
        assert len(x[1]) == 1
        assert len(y[1]) == 1
        assert isinstance(x[1], set)
        assert isinstance(y[1], set)
        assert list(x[1])[0][0] == list(y[1])[0][0]
        assert list(x[1])[0][1] == list(y[1])[0][2]
        assert list(x[1])[0][2] == list(y[1])[0][1]

    def test_findpath_flooding_nr2(self):
        seq = "AAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUAC"
        s1 = "....(((.(.((((....))))).)))((((...((((........))))...))))..."
        s2 = "...(((((((......(((((......)))))......))))...)))............"

        md = RNA.md()
        minh = 4

        crn = []
       
        lmp = get_fpath_flooding_cache(seq, s1, s2, md, minh)

        print()
        print(seq)
        print(s1)
        print(s2)
        print(' -- ')

        self._print_multi_path(seq, s1, s2, seq, s1, s2)
 
    def test_findpath_flooding_random(self):

        #seq = "AGCCCCGCGGUUACUCUCGAUUCCGUACUCAACUUUGUGCUUGCGUAUGAAUCUCCUUAACGGAUUCGAACCCCAAGCCGUUGCUCAAGCGCCCGCGCCUUAGCUUUUGGUUGUGUGUCCGCACAAAUGUAUCAGGCUUCGUCGACUUGG" 
        #ss1 = ".((...((((((....((((.(((((.......(((((((....)))))))........))))).))))......)))))).)).((((((.(...((((.(((.((((..(((......))))))).)).).))))...).)).))))."
        #ss2 = "(((...(((((.....((((.(((((...........(((....)))............))))).)))).......))))).)))((((((.....((((..((.(....((((((....))))))).))...)))).....)).))))."
        #seq = "GGGGUCGCAGCGCCCACGGAGACGGAUCAAUAAUACCAUUUGGGCCAUCACCGUCGCGAUUCCCCGUAAACCACCGAUUUCGGCUAAGAGUCCCUACAGUGACAUAAAACCUAACCGUAAUCGGAUUUCUUCCGUCAAUGUGGCACUAGC"
        #ss1 = "((((...(.(((...((((.((.((.((((.........)))).)).)).))))))))...))))((...((((......(((..((((((((.(((.((.............)).)))...)))).))))))).....)))).))...."
        #ss2 = "((((.....(((...((((.((.((.....(((......)))..)).)).)))))))....))))((...((((......(((..((((((((.(((.((.............)).)))...)))).))))))).....)))).))...."
        #seq = "CUACGAUGCCGCCACCAUCUGCACACGUGAGCCAAAGAACUAAGCUCUUUUUGCGUAGAAUUGGCAAAGCUUUAGGCCGGGUUCGUGAGACUCCUUUACGUGUAAUUGCUGCAGAAUAGUCCAAAACAGAAUUCCCCGCAUAUCAAUCGG"
        #ss1 = "((.....((.((((...(((((.((...((((...........))))....)).)))))..))))...))...))((.((((((((..((((........((((.....))))....))))....)).)))..))).))..........."
        #ss2 = ".(.....((.((((....((((.((...((((...........))))....)).))))...))))...)).....).((((((((((.((((..(((...(((.......)))))).))))))...).)))...))))............"
        #seq = "AUAUAUGGGGACCGUAGCAAGACUAUUCACGGCAUUCUUCUCUGCAAAUGGUUAACCAAUUUUGCAGUCUGGCUGACAGCUCGAUUUGCAUGCACAGCCCGACCCCGUUGGACGGGUGGGAUCUCACCUUAUCUCAGCCGCCAUAUUCUG"
        #ss1 = "....((((((.((((.............)))).........((((((.(((....)))...))))))((.(((((.((((.......)).))..))))).)))))))).((.(((.((((((........)))))).)))))........"
        #ss2 = "......((((...((((.....))))...............((((((((((....)))..)))))))((.(((((.((((.......)).))..))))).))))))...((.(((.((((((........)))))).)))))........"
        #seq = "CUUUUCUCUUAGUUAUGACUUUCCCAUCGAAACACGUGCAGUAUAGAAAUCCAAAUGCGCUUUGCGGGGCAGGGUACCGUGCCUUUAGCUCGAGCUGCUGUCCUCGAAAAAAGUCUAGUUUGGCGCGGAUCCGACCAUGGCAUCCGCUGG"
        #ss1 = "...((((............((((.....))))............))))..((((((..(((((.(((((((((((..((.((.....)).)).))).)))))).))...)))))...)))))).((((((((......)).))))))..."
        #ss2 = "......(((((.(((((..((((.....))))..)))).).)).)))...((((((..((((..(((((((((((..((.((.....)).)).))).))).)))))....))))...)))))).((((((((......)).))))))..."
        seq = "AGCUACGCCGUGUGCCCUAGGGCUUUCGUUGAUCCCACAGCCGCCAAACGACGGGUUAUACUGGGCAAUGAGAUCUCUAAUGGUGGAUACAUGAGGGAUUACUGGCGCGAGGAGUUGAAUUGAGUAGUGUAUCCUCAUUUCUUGUAUUGC"
        ss1 = ".....((((((...((((((((.((((((((..((((.....(((........))).....)))))))))))).))))).(.(((....))).))))...)).))))((((((((.((...((........)).))))))))))......"
        ss2 = ".....((((((......(((((.((((((((..((((.....(((........))).....)))))))))))).)))))))))))((((((((((((((((((....(((....))).....))))))...)))))).....)))))).."

        md = RNA.md()
        minh = 3

        crn = []
       

        print()
        print(seq)
        print(ss1)
        print(ss2)
        print(' -=- multi-path -=- ')

        _ = get_fpath_flooding_cache(seq, ss1, ss2, md, minh)

        self._print_multi_path(seq, ss1, ss2, seq, ss1, ss2)


    def _print_multi_path(self, seq, ss1, ss2, startseq, start, stop):
        lmp = get_fpath_flooding_cache(seq, ss1, ss2, None, None)
        if lmp[0] is None:
            #print(lmp[1])
            if isinstance(lmp[1], set):
                for (lsq, pm1, pm2) in lmp[1]: 
                    #print(lsq, pm1, pm2)
                    self._print_multi_path(lsq, pm1, pm2, startseq, start, stop)
                    start = apply_bp_change(startseq, start, stop, lsq, pm1, pm2)
            else:
                for (pm1, pm2) in lmp[1]: 
                    #print(pm1, pm2)
                    self._print_multi_path(seq, pm1, pm2, startseq, start, None)
                    start = apply_bp_change(startseq, start, None, seq, pm1, pm2)
        else:
            stop = apply_bp_change(startseq, start, stop, seq, ss1, ss2)
            revp = get_fpath_flooding_cache(seq, ss2, ss1, None, None)
            print(start, '\\/ fowd', lmp)
            print(stop, 'back /\\', revp)


if __name__ == '__main__':
    unittest.main()
