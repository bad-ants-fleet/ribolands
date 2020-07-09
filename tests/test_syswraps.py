#!/usr/bin/env python3

import os
import gzip
import unittest

import ribolands.syswraps as rsys

# Assumes these programs to be executable.
RNAsubopt = 'RNAsubopt'
barriers = 'barriers'
treekin = 'treekin'

def missing(program):
    return rsys.which(program) is None

@unittest.skip(f"need to implement io handling!")
class Test_Workflow(unittest.TestCase):
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
    Pipe = rsys.Workflow(seq, vrna_md, name = name)

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

class Test_pipeline(unittest.TestCase):
    def setUp(self):
        self.testname = 'ribolands_testsuite'
        self.files = set()
        # Check if files would be overwritten, if everything is fine, add them to
        # self.files. They get cleaned up after the tests have completed.
        for f in self._RiboTestFiles(self.testname):
            self.files.add(f)

    def tearDown(self):
        # remove files that might have been written to disc
        for fname in self.files:
            if os.path.exists(fname):
                os.remove(fname)

    def test_sys_subopt_range(self):
        seq = 'UAACUCACAAUGGUUGCAAA'

        with self.assertRaises(rsys.ExecError):
            rsys.sys_subopt_range(seq, RNAsubopt='x')

        # Testing defaults, adapt here if they change!
        assert rsys.sys_subopt_range(seq) == (30, 2937)
        assert rsys.sys_subopt_range(seq, nos = 1) == (0, 1)
        assert rsys.sys_subopt_range(seq, nos = 20) == (3.1, 22)
        # Testing hidden feature: return number of structs for given energy:
        assert rsys.sys_subopt_range(seq, nos = 0, maxe = 8.5) == (8.5, 390)

    def test_bugs_subopt_range(self):
        # These bugs have been fixed by increasing the starting energy range for
        # RNAsubopt computation from 2kcal/mol to 5kcal/mol.
        seq = 'UAACUCACAAUGGUUGCAAA'

        # Those results may be off because RNAsubopt returns a wrong density of states output.
        assert rsys.sys_subopt_range(seq, nos = 1, circ = True) == (0.0, 1)
        assert rsys.sys_subopt_range(seq, nos = 2, circ = True) == (4.6, 2)
        assert rsys.sys_subopt_range(seq, nos = 10, circ = True) == (6.6, 10)
        assert rsys.sys_subopt_range(seq, nos = 20, circ = True) == (8.4, 20)
        assert rsys.sys_subopt_range(seq, nos = 0, maxe = 8.5, circ = True) == (8.5, 20)
        assert rsys.sys_subopt_range(seq, nos = 20, circ = True, temp = 10) == (6.0, 20)

        seq = 'GGGUCGCCGUUACAUAGACCCUGCAACUAU'
        assert rsys.sys_subopt_range(seq, nos = 7100000, maxe = 20.00) == (20.0, 60872)

    def test_sys_suboptimals(self):
        seq = 'UAACUCACAAUGGUUGCAAA'
        name = self.testname
        ref_name = 'tests/files/' + self.testname + '.spt.gz'
        spt_name = rsys.sys_suboptimals(name, seq, ener = 1.0)

        # If this breaks, pass the tearDown function to keep temporary files.
        with gzip.open(spt_name) as spt, gzip.open(ref_name) as ref:
            for r in ref:
                s = spt.readline()
                self.assertEqual(s, r)

    @unittest.skipIf(missing(barriers), reason = f"Cannot find {barriers} executable!")
    def test_sys_barriers(self):
        seq = 'UAACUCACAAUGGUUGCAAA'
        name = self.testname
        sfile = 'tests/files/' + self.testname + '.spt.gz'
        ref_bfile = 'tests/files/' + self.testname + '.bar'
        ref_efile = 'tests/files/' + self.testname + '.err'
        ref_rfile = 'tests/files/' + self.testname + '.rts'

        [bofile, befile, brfile, bbfile, bpfile, bmfile] = rsys.sys_barriers_180(name, sfile)

        self.assertFileWeakEqual(bofile, ref_bfile)
        self.assertFileWeakEqual(brfile, ref_rfile)

    @unittest.skipIf(missing(treekin), reason = f"Cannot find {treekin} executable!")
    def test_sys_treekin(self):
        seq = 'UAACUCACAAUGGUUGCAAA'
        name = self.testname
        bofile = 'tests/files/' + self.testname + '.bar'
        brfile = 'tests/files/' + self.testname + '.rts'
        ref_tfile = 'tests/files/' + self.testname + '.tkn'

        with self.assertRaises(rsys.SubprocessError): # Not ergodic!
            tofile, tefile = rsys.sys_treekin_051(name, brfile, bofile = bofile)

    @unittest.skipIf(missing(treekin) or missing(barriers), 
            reason = f"Cannot find {treekin} or {barriers} executable!")
    def test_full_workflow(self):
        name = 'full_workflow_1'
        seq = 'UAACUCACAAUGGUUGCAAA'
        ref_tfile = 'tests/files/' + name + '.tkn'
        ref_efile = 'tests/files/' + name + '.err'

        for f in self._RiboTestFiles(name):
            self.files.add(f)

        sfile = rsys.sys_suboptimals(name, seq, ener = 5.0)
        [bofile, befile, brfile, bbfile, bpfile, bmfile] = rsys.sys_barriers_180(name, sfile, maxn = 7)
        tofile, tefile = rsys.sys_treekin_051(name, bbfile, bofile = bofile, t0 = 0, ti = 2, t8 = 1e10)

        # NOTE so many treekin changes...
        #self.assertFileWeakEqual(tfile, ref_tfile)

    def _RiboTestFiles(self, prefix):
        files = [
            prefix + '.spt.gz',
            prefix + '.bar',
            prefix + '.err',
            prefix + '_rates.txt',
            prefix + '_rates.bin',
            prefix + '.ps',
            prefix + '.tkn']
        for fname in files:
            if os.path.exists(fname):
                raise Exception(f"Attempting to overwrite temporary file: {fname}")
        return files

    def assertFileWeakEqual(self, test_file, ref_file):
        """Read two files and compare them line-by-line.

        Weak, beacuse it skipps lines staring with '%' or '#'.
        """
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
