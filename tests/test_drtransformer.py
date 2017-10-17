#!/usr/bin/env python

import sys
import unittest

import RNA
import ribolands.trafo as trafo

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
