#!/usr/bin/env python
#
# Written by Stefan Badelt (stef@tbi.univie.ac.at)
#
# ribolands.syswraps testing
#

import os
import gzip
import filecmp
import unittest
import ribolands.syswraps as rsys
from scipy import constants

class Test_complete_pipeline(unittest.TestCase):

  def setUp(self):
    self.testname = 'ribolands_testsuite' 
    self.files = set()
    # Check if files would be overwritten, if everything is fine, add them to
    # self.files. They get cleaned up after the tests have completed.
    for f in self._RiboTestFiles(self.testname):
      self.files.add(f)

  def tearDown(self):
    # remove files that might have been written to disc
    for fname in self.files :
      if os.path.exists(fname) : 
        os.remove(fname)

  def test_sys_subopt_range(self):
    seq  = 'UAACUCACAAUGGUUGCAAA'

    with self.assertRaises(rsys.ExecError):
      rsys.sys_subopt_range(seq, RNAsubopt='x')

    # Testing defaults, adapt here if they change!
    self.assertEqual(rsys.sys_subopt_range(seq), (30, 2937))

    self.assertEqual(rsys.sys_subopt_range(seq, nos=20), (3.1, 22))

    self.assertEqual(rsys.sys_subopt_range(seq, nos=1), (0, 1))

    # Testing hidden feature: return number of structs for given energy:
    self.assertEqual(rsys.sys_subopt_range(seq, nos=0, maxe=8.5), (8.5, 390))

  def test_bugs_subopt_range(self):
    # These bugs have been fixed by increasing the starting energy range for 
    # RNAsubopt computation from 2kcal/mol to 5kcal/mol.
    seq  = 'UAACUCACAAUGGUUGCAAA'

    # Returns energy-range 2.0 kcal/mol with 1 structure, because the secondary
    # structure space starts at around 4 kcal/mol.
    #self.assertEqual(rsys.sys_subopt_range(seq, nos=20, circ=True), (2.0, 1))
    self.assertEqual(rsys.sys_subopt_range(seq, nos=20, circ=True), (8.4, 20))

    # Hack the function to show that there exists a better result!
    self.assertEqual(rsys.sys_subopt_range(seq, nos=0, maxe=8.5, circ=True), (8.5, 20))

    # Interestingly, this effect gets reduced when reducing the temperature,
    # but it terminates at 10 structures, instead of moving on ...
    #self.assertEqual(rsys.sys_subopt_range(seq, nos=20, circ=True, temp=10), (4.7, 10))
    self.assertEqual(rsys.sys_subopt_range(seq, nos=20, circ=True, temp=10), (6.0, 20))

    seq2 = 'GGGUCGCCGUUACAUAGACCCUGCAACUAU'
    self.assertEqual(rsys.sys_subopt_range(seq2, nos=7100000, maxe=20.00), (20.0, 60872))

  def test_sys_suboptimals(self):
    seq  = 'UAACUCACAAUGGUUGCAAA'
    name = self.testname
    ref_name = 'tests/files/' + self.testname + '.spt.gz'
    spt_name = rsys.sys_suboptimals(name, seq, ener=1.0)

    # If this breaks, pass the tearDown function to keep temporary files.
    with gzip.open(spt_name) as spt, gzip.open(ref_name) as ref :
      for r in ref :
        s = spt.readline()
        self.assertEqual(s, r)

  def test_sys_barriers(self):
    seq  = 'UAACUCACAAUGGUUGCAAA'
    name = self.testname
    sfile = 'tests/files/' + self.testname + '.spt.gz'
    ref_bfile = 'tests/files/' + self.testname + '.bar'
    ref_efile = 'tests/files/' + self.testname + '.err'
    ref_rfile = 'tests/files/' + self.testname + '.rts'
    ref_psfile = 'tests/files/' + self.testname + '.ps'

    [sfile, bfile, efile, rfile, psfile] = rsys.sys_barriers(name, seq, sfile)

    self.assertFileWeakEqual(bfile, ref_bfile)
    self.assertFileWeakEqual(efile, ref_efile)
    self.assertFileWeakEqual(rfile, ref_rfile)
    self.assertFileWeakEqual(psfile, ref_psfile)

  def test_sys_treekin(self):
    seq  = 'UAACUCACAAUGGUUGCAAA'
    name = self.testname
    bfile = 'tests/files/' + self.testname + '.bar'
    rfile = 'tests/files/' + self.testname + '.rts'
    ref_tfile = 'tests/files/' + self.testname + '.tkn'

    with self.assertRaises(rsys.SubprocessError):
      tfile, efile = rsys.sys_treekin(name, seq, bfile, rfile)

  def test_full_workflow(self):
    name = 'full_workflow_1'
    seq  = 'UAACUCACAAUGGUUGCAAA'
    ref_tfile = 'tests/files/' + name + '.tkn'
    ref_efile = 'tests/files/' + name + '.err'

    for f in self._RiboTestFiles(name):
      self.files.add(f)

    sfile = rsys.sys_suboptimals(name, seq, ener=4.0)
    [sfile, bfile, efile, rfile, psfile] = rsys.sys_barriers(name, seq, sfile)
    tfile, efile = rsys.sys_treekin(name, seq, bfile, rfile, t0=1e-6, ti=2, t8=1e10)

    self.assertFileWeakEqual(tfile, ref_tfile)
    #self.assertFileWeakEqual(efile, ref_efile)

  def _RiboTestFiles(self, prefix):
    files = [
        prefix + '.spt.gz',
        prefix + '.bar',
        prefix + '.err',
        prefix + '.rts',
        prefix + '.ps',
        prefix + '.tkn' ]

    for fname in files :
      if os.path.exists(fname) : 
        raise Exception("attempting to overwrite temporary file", fname)

    return files

  def assertFileWeakEqual(self, test_file, ref_file):
    """Read two files and compare them line-by-line. 

    Weak, beacuse it skipps lines staring with '%' or '#'.
    """
    with open(test_file) as tst, open(ref_file) as ref :
      for r in ref :
        t = tst.readline()
        # A comment that may contain a current timestamp
        if t[0] == '%': continue
        if t[0] == '#': continue
        self.assertEqual(t, r)

