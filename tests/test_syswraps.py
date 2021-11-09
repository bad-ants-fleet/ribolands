#!/usr/bin/env python3

import os
import gzip
import unittest
from shutil import which

from ribolands.syswraps import (Workflow, ExecError, SubprocessError,
                                sys_treekin,
                                sys_barriers,
                                sys_kinfold, # untested!
                                sys_suboptimals,
                                sys_subopt_range)

# Assumes these programs to be executable.
RNAsubopt = 'RNAsubopt'
barriers = 'barriers'
treekin = 'treekin'

def missing(program):
    return which(program) is None

class WorkflowTest(unittest.TestCase):
    def test_workflow(self):
        import RNA 
        from ribolands.parser import parse_barriers

        vrna_md = RNA.md()
        vrna_md.noLP = 0
        vrna_md.temperature = 25
        vrna_md.dangles = 2
        vrna_md.logML = 0
        vrna_md.special_hp = 1
        vrna_md.noGU = 0

        seq = 'CUCGUCGCCUUAAUCCAGUGCGGGCGCUAGACAUCUAGUUAUCGCCGCAA'
        seq = 'AAAGCCGCCUUAAGCCUACUUAGAUGGAAGUGACGUACGGGUAUUGGUACACGAUUUUAC'
        name = 'test'
        Pipe = Workflow(seq, vrna_md, name = name)

        Pipe.force = True
        Pipe.verbose = True
        Pipe.outdir = 'tests/files'
        Pipe.RNAsubopt = 'RNAsubopt'
        Pipe.barriers = 'barriers'
        Pipe.treekin = 'treekin'
        Pipe.Kinfold = 'Kinfold'
        Pipe.zip_suboptimals = False

        Pipe.find_subopt_range(nos = 10000, maxe = 20)
        sofile = Pipe.call_RNAsubopt()
        bofile, befile, brfile, bbfile, *_ = Pipe.call_barriers(
                minh = 3, maxn = 20, plot = False, bsize = True,
                connected = True, force = True)
        tofile, _ = Pipe.call_treekin(p0 = ['8=1'], useplusI = True)

        lmins = parse_barriers(bofile, return_tuple = True) 
        assert lmins[0] == Pipe.sequence
        for fn in Pipe.files:
            if fn is None:
                continue
            os.remove(fn)

class Test_pipeline(unittest.TestCase):
    def setUp(self):
        self.testname = 'ribolands_testsuite'
        self.files = set()
        # Check if files would be overwritten, if everything is fine, add them to
        # self.files. They get cleaned up after the tests have completed.
        for fname in self._RiboTestFiles(self.testname):
            if os.path.exists(fname):
                raise AssertionError(f"Temporary file exists: {fname}")
            self.files.add(fname)

    def tearDown(self):
        # remove files that might have been written to disc
        for fname in self.files:
            if os.path.exists(fname):
                os.remove(fname)

    def _RiboTestFiles(self, prefix):
        return [prefix + '.spt.gz',
                prefix + '_barriers.bar',
                prefix + '_barriers.err',
                prefix + '_barriers_rates.txt',
                prefix + '_barriers_rates.bin',
                prefix + '_treekin.nxy',
                prefix + '_treekin.err']

    def test_sys_subopt_range(self):
        seq = 'UAACUCACAAUGGUUGCAAA'
        with self.assertRaises(ExecError):
            sys_subopt_range(seq, RNAsubopt = 'x')
        # Testing defaults, adapt here if they change!
        assert sys_subopt_range(seq) == (30, 2937)
        assert sys_subopt_range(seq, nos = 1) == (0, 1)
        assert sys_subopt_range(seq, nos = 20) == (3.1, 22)
        # Testing hidden feature: return number of structs for given energy:
        assert sys_subopt_range(seq, nos = 0, maxe = 8.5) == (8.5, 390)

    def test_bugs_subopt_range(self):
        # These bugs have been fixed by increasing the starting energy range for
        # RNAsubopt computation from 2kcal/mol to 5kcal/mol.
        seq = 'UAACUCACAAUGGUUGCAAA'

        # Those results may be off because RNAsubopt returns a wrong density of states output.
        assert sys_subopt_range(seq, nos = 1, circ = True) == (0.0, 1)
        assert sys_subopt_range(seq, nos = 2, circ = True) == (4.6, 2)
        assert sys_subopt_range(seq, nos = 10, circ = True) == (6.6, 10)
        assert sys_subopt_range(seq, nos = 20, circ = True) == (8.4, 20)
        assert sys_subopt_range(seq, nos = 0, maxe = 8.5, circ = True) == (8.5, 20)
        assert sys_subopt_range(seq, nos = 20, circ = True, temp = 10) == (6.0, 20)

        seq = 'GGGUCGCCGUUACAUAGACCCUGCAACUAU'
        assert sys_subopt_range(seq, nos = 7100000, maxe = 20.00) == (20.0, 60872)

    def test_sys_suboptimals(self):
        seq = 'UAACUCACAAUGGUUGCAAA'
        name = self.testname
        ref_name = 'tests/files/' + self.testname + '.spt.gz'
        spt_name = sys_suboptimals(name, seq, ener = 1.0)
        with gzip.open(spt_name) as spt, gzip.open(ref_name) as ref:
            for r in ref:
                s = spt.readline()
                self.assertEqual(s, r)
        os.remove(spt_name)

    @unittest.skipIf(missing(barriers), reason = f"Cannot find {barriers} executable!")
    def test_sys_barriers(self):
        seq = 'UAACUCACAAUGGUUGCAAA'
        name = self.testname
        sfile = 'tests/files/' + self.testname + '.spt.gz'
        ref_bfile = 'tests/files/' + self.testname + '_barriers.bar'
        ref_efile = 'tests/files/' + self.testname + '_barriers.err'
        ref_rfile = 'tests/files/' + self.testname + '_barriers_rates.txt'

        [bofile, befile, brfile, bbfile, bpfile, bmfile] = sys_barriers(name, sfile)

        self.assertFileWeakEqual(bofile, ref_bfile)
        self.assertFileWeakEqual(brfile, ref_rfile)
        os.remove(bofile)
        os.remove(befile)
        os.remove(brfile)
        os.remove(bbfile)
        assert bpfile is None
        assert bmfile is None

    @unittest.skipIf(missing(treekin), reason = f"Cannot find {treekin} executable!")
    def test_sys_treekin(self):
        seq = 'UAACUCACAAUGGUUGCAAA'
        name = self.testname
        bofile = 'tests/files/' + self.testname + '_barriers.bar'
        brfile = 'tests/files/' + self.testname + '_barriers_rates.txt'
        ref_tfile = 'tests/files/' + self.testname + '_treekin.nxy'
        with self.assertRaises(SubprocessError): # Not ergodic!
            tofile, tefile = sys_treekin(name, brfile, bofile = bofile)

    @unittest.skipIf(missing(treekin) or missing(barriers), 
            reason = f"Cannot find {treekin} or {barriers} executable!")
    def test_full_workflow(self):
        name = 'full_workflow'
        seq = 'UAACUCACAAUGGUUGCAAA'
        ref_tfile = 'tests/files/' + name + '_treekin.nxy'
        ref_efile = 'tests/files/' + name + '_treekin.err'

        for f in self._RiboTestFiles(name):
            self.files.add(f)

        sfile = sys_suboptimals(name, seq, ener = 5.0)
        [bofile, befile, brfile, bbfile, bpfile, bmfile] = sys_barriers(name, sfile, maxn = 7)
        tofile, tefile = sys_treekin(name, bbfile, bofile = bofile, t0 = 0, ti = 2, t8 = 1e10)

    def assertFileWeakEqual(self, test_file, ref_file):
        """ Read two files and compare them line-by-line.

        Weak, beacuse it skipps lines staring with '%' or '#'.
        """
        #print(f'comparing {test_file=} {ref_file=}')
        with open(test_file) as tst, open(ref_file) as ref:
            for r in ref:
                t = tst.readline()
                # A comment that may contain a current timestamp
                if t[0] == '%':
                    continue
                if t[0] == '#':
                    continue
                self.assertEqual(t, r)

if __name__ == '__main__':
    unittest.main()
