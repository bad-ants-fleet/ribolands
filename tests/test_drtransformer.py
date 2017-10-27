#!/usr/bin/env python

import sys
import unittest

import RNA
import ribolands.trafo as trafo
from ribolands.syswraps import sys_treekin

class Test_ConformationGraph(unittest.TestCase):
  def setUp(self):
    pass

  def dont_test_minitrafo(self):
    # remove the dont_ for testing, but beware it writes files ...
    fullseq = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCA"
    vrna_md = RNA.md()

    CG = trafo.ConformationGraph(fullseq, vrna_md)

    self.assertEqual(CG.nodes(), [])

    CG._transcript_length = 1
    CG.expand()
    self.assertEqual(len(CG), 1)
    self.assertEqual(CG.nodes(), ['.' * len(fullseq)])
    self.assertEqual(CG._nodeid, 1)

    [bfile, rfile, p0, nlist] = CG.get_simulation_files_tkn('rudi')

    stepsize = 2

    for i in range(2, len(fullseq), stepsize):
      seq = fullseq[0:i]
      CG.expand(extend=stepsize)
      #if i == 10:
      #  CG.logfile = sys.stdout
      [bfile, rfile, p0, nlist] = CG.get_simulation_files_tkn('rudi')
      if len(nlist) == 1 :
        CG._total_time += 0.2
      else :
        bfile = None # sometimes bfile causes a segfault, so let's leave it out.
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

    CG = trafo.ConformationGraph(fullseq, vrna_md)

    self.assertEqual(CG.nodes(), [])

    CG.expand()
    self.assertEqual(len(CG), 1)
    self.assertEqual(CG.nodes(), ['.' * len(fullseq)])
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
        [".(((......(((.((((((.....))))))..)))......))).....", 0.25, True] ]
    for (ss,occ,active) in ess :
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
        [".(((......(((.((((((.....))))))..)))......))).....", 0.25, True] ]
    for (ss,occ,active) in ess :
      self.assertTrue(CG.has_node(ss))
      self.assertEqual(CG.node[ss]['occupancy'], occ)
      self.assertEqual(CG.node[ss]['active'], active)

    #CG.get_simulation_files_tkn('expand')

    CG.coarse_grain(minh=4.3)
    ess = [
        ["..........((((((..((((...((....))...)))).))))))...", 0.25, True], 
        ["........................((((((........))))))......", 0.00, True],
        [".......(((((((...)))))))((((((........))))))......", 0.50, True],
        [".........((((((....).))))).(((.......)))..........", 0.00, True],
        ["........................((((((.......))).)))......", 0.00, False],
        [".......(((((((...)))))))((((((.......))).)))......", 0.00, False],
        [".....(((..(((.((((((.....))))))..))).....)))......", 0.25, True],
        [".(((......(((.((((((.....))))))..)))......))).....", 0.00, False] ]
    for (ss,occ,active) in ess :
      self.assertTrue(CG.has_node(ss))
      self.assertEqual(CG.node[ss]['occupancy'], occ)
      self.assertEqual(CG.node[ss]['active'], active)

    CG.prune()
    #CG.logfile = sys.stdout
    #CG.get_simulation_files_tkn('expand_again')
    ess = [
        ["..........((((((..((((...((....))...)))).))))))...", 0.25, True], 
        ["........................((((((........))))))......", 0.00, False],
        [".......(((((((...)))))))((((((........))))))......", 0.50, True],
        [".........((((((....).))))).(((.......)))..........", 0.00, False],
        ["........................((((((.......))).)))......", 0.00, False],
        [".......(((((((...)))))))((((((.......))).)))......", 0.00, False],
        [".....(((..(((.((((((.....))))))..))).....)))......", 0.25, True],
        [".(((......(((.((((((.....))))))..)))......))).....", 0.00, False] ]
    for (ss,occ,active) in ess :
      self.assertTrue(CG.has_node(ss))
      self.assertEqual(CG.node[ss]['occupancy'], occ)
      self.assertEqual(CG.node[ss]['active'], active)
   
  def initialize_CG(self, seq, sss):
    fullseq = seq
    vrna_md = RNA.md()

    CG = trafo.ConformationGraph(fullseq, vrna_md)
    CG._transcript_length = len(seq)

    for e, s1 in enumerate(sss, 1) :
      if not CG.has_node(s1):
        en = round(CG._fold_compound.eval_structure(s1), 2)
        CG.add_node(s1, energy=en, occupancy=1.0/len(sss), 
            identity=CG._nodeid, active=True, last_seen=0)
        CG._nodeid += 1
      for s2 in sss[e:]:
        if not CG.has_node(s2):
          en = round(CG._fold_compound.eval_structure(s2), 2)
          CG.add_node(s2, energy=en, occupancy=1.0/len(sss), 
              identity=CG._nodeid, active=True, last_seen=0)
          CG._nodeid += 1
        assert CG.add_transition_edge(s1, s2)

    return CG



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

    self.assertEqual(nbr, '.....(((((......)).)))((((((((....))))....))))..')
    self.assertEqual(ext,                'CUCGUCGNNNCGGGNNNCCGC')
    self.assertEqual(ext_moves[ext][1],  '.....((...))((...))..')
    self.assertEqual(ext_moves[ext][0],  set())

  def test_open_breathing_helices(self):
    se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGC"
    ss = ".....(((((......)).)))((((((((....))))....)))).."

    out = trafo.open_breathing_helices(se, ss, free = 6)

    res = [
      "........((......))....((((((((....))))....))))..",
      ".....(((((......)).)))....((((....))))..........",
      "........((......))........((((....)))).........."]

    self.assertEqual(sorted(out), sorted(res))

    out = trafo.open_breathing_helices(se, ss, free = 8)

    res = [
      "......................((((((((....))))....))))..",
      ".....(((((......)).)))....((((....))))..........",
      "..........................((((....)))).........."]

    self.assertEqual(sorted(out), sorted(res))

  def test_open_breathing_helices_multi(self):
    se = "CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCG"
    ss = "..((.(((((......)).)))((((((((....))))....)))).))"

    out = trafo.open_breathing_helices(se, ss, free = 6)
    res = [".....(((((......)).)))((((((((....))))....))))..."]
    self.assertEqual(sorted(out), sorted(res))

    out = trafo.open_breathing_helices(se, ss, free = 7)
    res = [
        "........((......))....((((((((....))))....))))...",
        ".....(((((......)).)))....((((....))))...........",
        "........((......))........((((....))))..........."]
    self.assertEqual(sorted(out), sorted(res))



