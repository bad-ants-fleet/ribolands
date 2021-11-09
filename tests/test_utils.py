#!/usr/bin/env python

import unittest
import subprocess as sub

from ribolands.utils import RiboUtilsError
import ribolands.utils as rutils

class test_io(unittest.TestCase):
    # fasta-like input parser that takes single name/sequence input
    def test_vienna_stdin(self):
        # single sequence input
        regular_input = ['>testseq with additional info in header',
                         'CUUGUUCGAUGUUGUCUUUCACGAAUCCGGUCCGCAACCCAUUAUAAGCAGAUGUGUAG']

        fasta_input = ['>testseq with additional info in header',
                       'CUUGUUCGAUGU ',
                       'UGUCUUUCACGAAUCCGGUCC',
                       ' GCAACCCAUUAUAAGCAGAUGUGUAG',
                       '']
        result = (
            'testseq',
            'CUUGUUCGAUGUUGUCUUUCACGAAUCCGGUCCGCAACCCAUUAUAAGCAGAUGUGUAG')
        self.assertEqual(rutils.parse_vienna_stdin(regular_input), result)
        self.assertEqual(rutils.parse_vienna_stdin(fasta_input), result)

    def test_vienna_stdin_cofold(self):
        # cofold input with multiple cutpoints
        multifold_input = ['>testseq with additional info in header',
                           'CUUGUUCGAUGU& ',
                           'UGUCUUUCACGAAUCCGGUCC&GCAACCCAUUAUAAGCAGAUGUGUAG']

        result = (
            'testseq',
            'CUUGUUCGAUGU&UGUCUUUCACGAAUCCGGUCC&GCAACCCAUUAUAAGCAGAUGUGUAG')
        self.assertEqual(rutils.parse_vienna_stdin(multifold_input), result)

    def test_parse_ratefile(self):
        # TODO
        pass

class dot_bracket_modifications(unittest.TestCase):
    # Make a pair table representation of secondary structures.
    def test_pair_table(self):
        ss = '(((...)).((...)).((...)))'
        pt_b0 = [24, 7, 6, -1, -1, -1, 2, 1, -1, 15, 14, -1, -1, -1, 10, 9, -1, 23, 22, -1, -1, -1, 18, 17, 0]
        pt_b1 = [25, 25, 8, 7, 0, 0, 0, 3, 2, 0, 16, 15, 0, 0, 0, 11, 10, 0, 24, 23, 0, 0, 0, 19, 18, 1]
        assert pt_b0 == rutils.make_pair_table(ss)
        assert pt_b1 == rutils.make_pair_table(ss, base = 1)

    def test_pair_table_errors(self):
        ss = '(((...)).((&)).((...)))'
        with self.assertRaises(RiboUtilsError):
            rutils.make_pair_table(ss)
        pt_b0 = rutils.make_pair_table(ss, chars = list('&.'))
        pt_b1 = rutils.make_pair_table(ss, chars = list('&.'), base = 1)
        assert pt_b0 == [22, 7, 6, -1, -1, -1, 2, 1, -1, 13, 12, -1, 10, 9, -1, 21, 20, -1, -1, -1, 16, 15, 0]
        assert pt_b1 == [23, 23, 8, 7, 0, 0, 0, 3, 2, 0, 14, 13, 0, 11, 10, 0, 22, 21, 0, 0, 0, 17, 16, 1]

        ss = '(((...)).((.).((...)))'
        with self.assertRaises(RiboUtilsError):
            rutils.make_pair_table(ss)

        ss = '((...)).((.)).((...)))'
        with self.assertRaises(RiboUtilsError):
            rutils.make_pair_table(ss)

    def test_loop_index(self):
        ss = '.(((...)).((...))..(.(...)))'
        li = '0123333321455555411667777761'
        assert rutils.make_loop_index(ss) == list(map(int, li))
 
        ss = '.(((...)(...).((.(...))).)).'
        li = '0123333344444256677777652210'
        assert rutils.make_loop_index(ss) == list(map(int, li))
 
if __name__ == '__main__':
    unittest.main()
