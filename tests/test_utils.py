#!/usr/bin/env python

import sys
import unittest
import subprocess as sub
import ribolands.utils as ru

class Test_vienna_stdin(unittest.TestCase):
  """Standard fasta-like input parser that takes single name/sequence input."""

  def setUp(self):
    self.regular_input = [
        '>testseq with additional info in header',
        'CUUGUUCGAUGUUGUCUUUCACGAAUCCGGUCCGCAACCCAUUAUAAGCAGAUGUGUAG']
    self.weird_input  = [
        '>testseq with additional info in header',
        'CUUGUUCGAUGU ',
        'UGUCUUUCACGAAUCCGGUCC',
        ' GCAACCCAUUAUAAGCAGAUGUGUAG', 
        '']
    self.cofold_input = [
        '>testseq with additional info in header',
        'CUUGUUCGAUGU& ',
        'UGUCUUUCACGAAUCCGGUCC&GCAACCCAUUAUAAGCAGAUGUGUAG']

  def test_vienna_stdin(self):
    E_result = ('testseq', 'CUUGUUCGAUGUUGUCUUUCACGAAUCCGGUCCGCAACCCAUUAUAAGCAGAUGUGUAG')
    E_result_cofold = ('testseq', 'CUUGUUCGAUGU&UGUCUUUCACGAAUCCGGUCC&GCAACCCAUUAUAAGCAGAUGUGUAG')
    C_result = ru.parse_vienna_stdin(self.regular_input)
    C_result2 = ru.parse_vienna_stdin(self.weird_input)
    C_result_cofold = ru.parse_vienna_stdin(self.cofold_input)

    self.assertEqual(E_result, C_result, 'Vienna-type fasta parser')
    self.assertEqual(E_result, C_result2, 'Vienna-type fasta parser')
    self.assertEqual(E_result_cofold, C_result_cofold, 'Vienna-type fasta parser')

class Test_make_pair_table(unittest.TestCase):
  """Make a pair table representation of secondary structures."""

  def setUp(self):
    self.s1 = '(((...)).((...)).((...)))'
    self.s2 = '(((...)).((&)).((...)))'

  def test_pair_table(self):
    Exp_pt1_s0 = [24, 7, 6, -1, -1, -1, 2, 1, -1, 15, 14, -1, -1, -1, 10, 9, -1, 23, 22, -1, -1, -1, 18, 17, 0]
    Exp_pt1_s1 = [25, 25, 8, 7, 0, 0, 0, 3, 2, 0, 16, 15, 0, 0, 0, 11, 10, 0, 24, 23, 0, 0, 0, 19, 18, 1]
    Cmp_pt1_s0 = ru.make_pair_table(self.s1)
    Cmp_pt1_s1 = ru.make_pair_table(self.s1, base=1)

    self.assertEqual(Cmp_pt1_s0, Exp_pt1_s0, '0 based pairtable')
    self.assertEqual(Cmp_pt1_s1, Exp_pt1_s1, '1 based pairtable (C-style)')

if __name__ == '__main__':
  unittest.main()

