#!/usr/bin/env python

import sys
import unittest

import RNA
import ribolands.trafo as trafo
from ribolands.syswraps import sys_treekin

class Test_ConformationGraph(unittest.TestCase):

  def test_minitrafo(self):
    fullseq = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCA"
    vrna_md = RNA.md()

    CG = trafo.ConformationGraph(fullseq, vrna_md)

    self.assertEqual(CG.nodes(), [])

    CG._transcript_length = 1
    CG.expand()
    self.assertEqual(len(CG), 1)
    self.assertEqual(CG.nodes(), ['.' * len(fullseq)])
    self.assertEqual(CG._nodeid, 1)

    CG.get_simulation_files_tkn(None)
    [bfile, rfile, p0, nlist] = CG.get_simulation_files_tkn('rudi')

    stepsize = 2

    for i in range(2, len(fullseq), stepsize):
      seq = fullseq[0:i]
      CG.expand(extend=stepsize)
      if i == 10:
        CG.logfile = sys.stdout
      [bfile, rfile, p0, nlist] = CG.get_simulation_files_tkn('rudi')
      if len(nlist) == 1 :
        CG._total_time += 0.2
      else :
        bfile = None # sometimes bfile causes a segfault, so let's leave it out.
        tfile, _ = sys_treekin('rudi', seq, bfile, rfile, binrates=True,
            treekin='treekin', p0=p0, t0=0, ti=1.5, t8=0.2, mpack=False,
            exponent=False, useplusI=False, force=True, verb=False)

        time_inc, iterations = CG.update_occupancies_tkn(tfile, nlist)
        CG._total_time += time_inc

      dn,sr = CG.prune(nlist)

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

    print 'n', nbr
    print 'e', ext

  def test_conformationgraph(self):
    fullseq = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
    vrna_md = RNA.md()

    CG = trafo.ConformationGraph(fullseq, vrna_md)

    self.assertEqual(CG.nodes(), [])

