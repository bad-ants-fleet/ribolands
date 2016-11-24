#!/usr/bin/env python

import sys
import unittest
import subprocess as sub

import ribolands.utils as rutils

class vienna_stdin_test(unittest.TestCase):
  #fasta-like input parser that takes single name/sequence input
  def setUp(self):
    pass

  def test_vienna_stdin(self):
    # single sequence input
    regular_input = ['>testseq with additional info in header',
                     'CUUGUUCGAUGUUGUCUUUCACGAAUCCGGUCCGCAACCCAUUAUAAGCAGAUGUGUAG']

    fasta_input  = ['>testseq with additional info in header',
                    'CUUGUUCGAUGU ',
                    'UGUCUUUCACGAAUCCGGUCC',
                    ' GCAACCCAUUAUAAGCAGAUGUGUAG', 
                    '']
    result = ('testseq', 'CUUGUUCGAUGUUGUCUUUCACGAAUCCGGUCCGCAACCCAUUAUAAGCAGAUGUGUAG')
    self.assertEqual(rutils.parse_vienna_stdin(regular_input), result)
    self.assertEqual(rutils.parse_vienna_stdin(fasta_input), result)

  def test_fasta_stdin_(self):
    # cofold input with multiple cutpoints
    multifold_input = ['>testseq with additional info in header',
                    'CUUGUUCGAUGU& ',
                    'UGUCUUUCACGAAUCCGGUCC&GCAACCCAUUAUAAGCAGAUGUGUAG']

    result = ('testseq', 'CUUGUUCGAUGU&UGUCUUUCACGAAUCCGGUCC&GCAACCCAUUAUAAGCAGAUGUGUAG')
    self.assertEqual(rutils.parse_vienna_stdin(multifold_input), result)

class make_pair_table_test(unittest.TestCase):
  #Make a pair table representation of secondary structures.
  def setUp(self):
    pass

  def test_pair_table(self):
    s1 = '(((...)).((...)).((...)))'
    Exp_pt1_s0 = [24, 7, 6, -1, -1, -1, 2, 1, -1, 15, 14, -1, -1, -1, 10, 9, -1, 23, 22, -1, -1, -1, 18, 17, 0]
    Exp_pt1_s1 = [25, 25, 8, 7, 0, 0, 0, 3, 2, 0, 16, 15, 0, 0, 0, 11, 10, 0, 24, 23, 0, 0, 0, 19, 18, 1]
    Cmp_pt1_s0 = rutils.make_pair_table(s1)
    Cmp_pt1_s1 = rutils.make_pair_table(s1, base=1)

    self.assertEqual(Cmp_pt1_s0, Exp_pt1_s0, '0 based pairtable')
    self.assertEqual(Cmp_pt1_s1, Exp_pt1_s1, '1 based pairtable (C-style)')

  def test_cofold_pair_table(self):
    s2 = '(((...)).((&)).((...)))'
    pass

if __name__ == '__main__':
  unittest.main()

