#!/usr/bin/env python
#
#  Coded by: Stefan Badelt <stef@tbi.univie.ac.at>
#  University of Vienna, Department of Theoretical Chemistry
#
#  -*- Content -*-
#  *) systemcalls of RNAsubopt, barriers and treekin
#  *) most likely requires linux
#
#  -*- Style -*-
#  Use double quotes or '#' for comments, such that single quotes are available
#  for uncommenting large parts during testing
#
import os
import re
import sys
import gzip
import math
import subprocess as sub
from struct import pack, unpack, calcsize

def check_version(program, rv):
    if which(program) is None:
        if 'treekin' in program:
            raise ExecError(program,
                            "treekin", 'http://www.tbi.univie.ac.at/RNA/Treekin')
        elif 'barries' in program:
            raise ExecError(barriers, "barriers",
                            'http://www.tbi.univie.ac.at/RNA/Barriers')
        elif 'RNAsubopt' in program:
            raise ExecError(RNAsubopt, "RNAsubopt",
                            'http://www.tbi.univie.ac.at/RNA')
        else:
            raise ExecError(program)
    p, pv = sub.check_output([program, '--version']).split()
    pv = pv.decode()

    def versiontuple(rv):
        return tuple(map(int, (rv.split("."))))

    if versiontuple(pv) < versiontuple(rv):
        raise VersionError(program, pv, rv)


class SubprocessError(Exception):
    """Raise Error: Commandline call failed."""

    def __init__(self, code, call=None):
        self.message = "Process terminated with return code: {}.\n".format(
            code)
        if call:
            self.message += "|==== \n"
            self.message += "{}\n".format(call)
            self.message += "|====\n"

        super(SubprocessError, self).__init__(self.message)


class ExecError(Exception):
    """Raise Error: Executable not found."""

    def __init__(self, var, program=None, download=None):
        self.message = "{} is not executable.\n".format(var)
        self.message += "|==== \n"
        if program:
            self.message += "|You need to install *{}*.\n".format(program)
            if download:
                self.message += "|Download: {}\n".format(download)

        self.message += "|If the program is installed, "
        self.message += "make sure that the path you specified is executable.\n"
        self.message += "|====\n"

        super(ExecError, self).__init__(self.message)


class VersionError(Exception):
    """Raise Error: Update program to latest version."""

    def __init__(self, var, cv, mv, download=None):
        self.message = "{} has version {}.\n".format(var, cv)
        self.message += "|==== \n"
        self.message += "|You need to install version v{} or higher.\n".format(
            mv)
        if download:
            self.message += "|Download: {}\n".format(download)
        self.message += "|====\n"

        super(VersionError, self).__init__(self.message)


# **************** #
# public functions #
# ................ #

def sys_treekin(name, seq, bfile, rfile,
                treekin='treekin',
                p0=['1=0.5', '2=0.5'],
                t0=0,
                ti=1.2,
                t8=1e6,
                binrates=False,
                useplusI=False,
                exponent=False,
                mpack_method='',
                force=False,
                verb=False):
    """Perform a system-call of the program ``treekin``.

    Prints the results into files and returns the respective filenames. This
    wrapper will produce two files, one from ``STDIN`` and ``STDERR`` output.

    Args:
      name (str): Name of the sequence used to name the output file.
      seq (str): Nucleic acid sequence.
      bfile (str): Input filename for ``treekin``. Fileformat as produced by ``barriers``.
      rfile (str): Input rate matrix for ``treekin``. Fileformat as produced by ``barriers``.
      p0 ([str,...], optional): A list to specify an initial occupancy vector. Defaults to:
        ['1=0.5', '2=0.5']
      t0 (float, optional): Start time of simulation.
      ti (float, optional): Time increment of the treekin solver.
      t8 (float, optional): Stop time for the solver.
      binrates (bool, optional): Ratefile is in binary format. Defaults to False.
      useplusI (bool, optional): Use treekin method --useplusI. Defaults to False.
      exponent (bool, optional): Use treekin method --exponent. Defaults to False.
      mpack_method (str, optional): Set treekin parameter --mpack-method. Defaults to ''.
      force (bool, optional): Overwrite existing files. Defaults to False.
      verb (bool, optional): Print verbose information. Defaults to False.

    Raises:
      ExecError: Program does not exist.
      SuboprocessError: Program terminated with exit code: ...

    Returns:
      [str, str]: The name of the file containing STDOUT and STDERR of ``treekin`` call.
    """

    if which(treekin) is None:
        raise ExecError(treekin,
                        "treekin", 'http://www.tbi.univie.ac.at/RNA/Treekin')

    reg_flt = re.compile(b'[-+]?[0-9]*.?[0-9]+([eE][-+]?[0-9]+)?.')
    # http://www.regular-expressions.info/floatingpoint.html

    tfile = name + '.tkn'
    efile = name + '.err'

    if not force and os.path.exists(tfile):
        if verb:
            print("# {:s} <= Files exist".format(tfile))
        return tfile, efile

    # Unfortunately, running treekin with a single state leads to an error that
    # is printed to STDOUT instead of STDERR. The program, exits with success.
    # That's why we catch it first:
    if binrates:
        with open(rfile, 'rb') as rf:
            i, = unpack('i', rf.read(calcsize('i')))  # unsigned long, little-endian
        i -= 1
    else:
        with open(rfile) as rf:
            for i, _ in enumerate(rf):
                pass
    if i == 0:
        raise SubprocessError(None, 'No transition rates found.')

    treecall = [treekin, '--method', 'I']
    if useplusI:
        treecall.extend(['--useplusI'])
    if mpack_method:
        treecall.extend(['--mpack-method=' + mpack_method])
    treecall.extend(['--tinc', str(ti)])
    treecall.extend(['--t0', str(t0)])
    treecall.extend(('--t8', str(t8)))
    for p in p0:
        treecall.extend(('--p0', p))
    if bfile:
        treecall.extend(('--bar', bfile))
    if binrates:
        treecall.extend(['--bin'])
    if exponent:
        treecall.extend(['--exponent'])

    if verb:
        if bfile:
            print("# {} < {} 2> {} > {}".format(
                ' '.join(treecall), bfile, efile, tfile))
        else:
            print("# echo "" | {} 2> {} > {}".format(
                ' '.join(treecall), efile, tfile))

    # Do the simulation (catch treekin errors)
    with open(rfile, 'r') as rts, \
            open(tfile, 'w') as tkn, \
            open(efile, 'w') as err:
        proc = sub.Popen(treecall, stdin=rts, stdout=tkn, stderr=err)
        proc.communicate(None)
        if proc.returncode:
            call = "{} < {} 2> {} > {}".format(
                ' '.join(treecall), bfile, efile, tfile)
            raise SubprocessError(proc.returncode, call)

    # Adapt here to return exact simulation time and number of iterations
    if exponent or verb:
        lastlines = sub.check_output(['tail', '-2', tfile]).strip().split(b'\n')
        if not reg_flt.match(lastlines[0]):
            raise SubprocessError(None, "No output from treekin simulation")
        else:
            if reg_flt.match(lastlines[1]):
                time = float(lastlines[1].split()[0])
                iterations = 0
            else:
                time = float(lastlines[0].split()[0])
                iterations = int(lastlines[-1].split()[-1])
            if (abs(float(time) - t8) > t8 / 1000):
                raise SubprocessError(None, "Treekin terminated at the wrong time.")
            if verb:
                print("# Treekin stopped after {:d} iterations at time {:f}".format(
                    iterations, time))

    return tfile, efile


def sys_barriers(name, seq, sfile,
                 barriers='barriers',
                 minh=0.001,
                 maxn=50,
                 k0=1.0,
                 temp=37.0,
                 noLP=False,
                 moves='single-base-pair',
                 gzip=True,
                 rates=True,
                 binrates=False,
                 bsize=False,
                 circ=False,
                 saddle=False,
                 mfile='',
                 force=False,
                 verb=False):
    """Perform a system-call of the program ``barriers``.

    The print the results into a file and return the respective filename. This
    wrapper will return the output files, including ``STDIN`` and ``STDERR``.

    Args:
      name (str): Name of the sequence used to name the output file.
      seq (str): Nucleic acid sequence.
      sfile (str): Input filename for ``barriers`` as produced by ``RNAsubopt``.
      temp (float, optional): Specify temperature in Celsius.
      force (bool, optional): Overwrite existing files with the same name.
      verb (bool, optional): Print the current system-call to ``stderr``.

    Raises:
      ExecError: Program does not exist.
      SuboprocessError: Program terminated with exit code: ...

    Returns:
      [sfile, bfile, rfile, efile, pfile]: A list of produced files containing
        ``barriers`` results.
    """

    if which(barriers) is None:
        raise ExecError(barriers, "barriers",
                        'http://www.tbi.univie.ac.at/RNA/Barriers')

    if not sfile or not os.path.exists(sfile):
        sfile = sys_suboptimals(name, seq,
                                verb=verb, noLP=noLP, circ=circ, temp=temp, force=force)

    bfile = name + '.bar'
    efile = name + '.err'
    rfile = name + '.rts'
    msfile = name + '.ms'
    psfile = name + '.ps'

    if not force and \
            os.path.exists(bfile) and \
            os.path.exists(rfile) and \
            os.path.exists(psfile) and \
            os.path.exists(efile):
        if verb:
            print("#", bfile, "<= Files exist")
        return [sfile, bfile, efile, rfile, psfile, msfile]

    if gzip:
        barcall = ['zcat', sfile, '|', barriers]
    else:
        barcall = ['cat', sfile, '|', barriers]

    if noLP:
        barcall.extend(['-G', 'RNA-noLP'])
    else:
        barcall.extend(['-G', 'RNA'])

    if moves == 'single-base-pair':
        pass
    elif moves == 'shift':
        barcall.extend(['--moves=Shift'])
    else:
        raise ValueError("Invalid move-set in sys_barriers()")

    # buggy barriers
    if subopt_reaches_minh(sfile, minh):
        barcall.extend(["--minh", str(minh)])
        barcall.extend(["--max", str(int(maxn))])
    barcall.extend(["-T", str(temp)])

    if rates:
        barcall.extend(['--rates'])
    if bsize:
        barcall.extend(['--bsize'])
    if saddle:
        barcall.extend(['--saddle'])

    if mfile:
        barcall.extend(["--mapstruc", mfile])

    if verb:
        print('#', ' '.join(barcall), '2>', efile, '>', bfile)

    with open(sfile, 'r') as shandle, \
            open(bfile, 'w') as bhandle, \
            open(efile, 'w') as ehandle:
        proc = sub.Popen(' '.join(barcall),
                         stdin=shandle, stdout=bhandle, stderr=ehandle, shell=True)
        proc.communicate(None)
        if proc.returncode:
            call = "{} 2> {} > {}".format(
                ' '.join(barcall), efile, bfile)
            raise SubprocessError(proc.returncode, call)

    if rates:
        if k0 != 1.0:
            if binrates:
                with open('rates.bin', 'rb') as rf, open(rfile, 'wb') as new:
                    dim, = unpack('i', rf.read(4))
                    new.write(pack("i", dim))
                    for e in range(dim * dim):
                        r, = unpack('d', rf.read(8))
                        new.write(pack("d", r * k0))
                os.rename('rates.out', rfile + '2')
            else:
                with open('rates.bin', 'rb') as rf, open(rfile, 'w') as new:
                    dim, = unpack('i', rf.read(4))
                    rm = []
                    for e in range(dim):
                        col = []
                        for e in range(dim):
                            r, = unpack('d', rf.read(8))
                            col.append(r * k0)
                        rm.append(col)
                    for line in zip(*rm):
                        newline = "".join(map("{:10.4g}".format, line))
                        new.write(newline + "\n")
                os.remove('rates.out')
            os.remove('rates.bin')
        else:
            if binrates:
                os.rename('rates.bin', rfile)
                os.rename('rates.out', rfile + '2')
            else:
                os.remove('rates.bin')
                os.rename('rates.out', rfile)
        os.remove('treeR.ps')

    if mfile:
        os.rename('mapstruc.out', msfile)
    os.rename('tree.ps', psfile)
    #print(sfile, bfile, efile, rfile, psfile, msfile)

    return [sfile, bfile, efile, rfile, psfile, msfile]


def sys_suboptimals(name, seq,
                    RNAsubopt='RNAsubopt',
                    ener=None,
                    temp=37.0,
                    verb=False,
                    noLP=False,
                    circ=False,
                    opts=[],
                    sort=['|', 'sort', '-T', '/tmp', '-k3r', '-k2n'],
                    gzip=True,
                    force=False):
    """ **Perform a system-call of the program ``RNAsubopt``.**

    The print the results into a file and return the filename.

    :param name: Name of the sequence used to name the output file.
    :param seq: Nucleic acid sequence.
    :param RNAsubopt: path to executable
    :param ener: Specify energy range
    :param temp: Specify temperature in Celsius.
    :param noLP: exclude lonely-base-pairs in suboptimal structures
    :param circ: compute density of states
    :param opts: more options for ``RNAsubopt``

    :param force: Overwrite existing files with the same name.
    :param verb: Print verbose information as a comment '#'.

    :type name: string
    :type seq: string
    :type RNAsubopt: string
    :type ener: float
    :type temp: float
    :type noLP: bool
    :type circ: bool
    :type opts: list
    :type force: bool
    :type verb: bool

    :return: Fliename of the file containing ``RNAsubopt`` results
    :rtype: string
    """

    if which(RNAsubopt) is None:
        raise ExecError(RNAsubopt, "RNAsubopt",
                        'http://www.tbi.univie.ac.at/RNA')

    if ener is None:
        ener, nos = sys_subopt_range(seq, verb=verb,
                                     RNAsubopt=RNAsubopt, noLP=noLP, circ=circ, temp=temp)
        if verb:
            print("# Energy-Update: {:.2f} kcal/mol to compute {} sequences".format(
                ener, nos))

    sfile = name + '.spt.gz' if gzip else name + '.spt'

    if os.path.exists(sfile) and not force:
        if verb:
            print("#", sfile, "<= File exists")
        return sfile

    sptcall = [RNAsubopt, "-e {:.2f} -T {:.2f}".format(ener, temp)]
    if noLP:
        sptcall.append("--noLP")
    if circ:
        sptcall.append("--circ")
    sptcall.extend(opts)
    sptcall.extend(sort)

    if gzip:
        sptcall.extend(('|', 'gzip', '--best'))

    if verb:
        print("#", "echo \"{}\" | {} > {}".format(seq, ' '.join(sptcall), sfile))

    with open(sfile, 'w') as shandle:
        proc = sub.Popen([' '.join(sptcall)],
                         stdin=sub.PIPE, stdout=shandle, shell=True)
        proc.communicate(seq.encode())
        if proc.returncode:
            call = "echo \"{}\" | {} > {}".format(
                seq, ' '.join(sptcall), sfile)
            raise SubprocessError(proc.returncode, call)

    return sfile


def sys_subopt_range(seq,
                     RNAsubopt='RNAsubopt',
                     nos=5100000,
                     maxe=30.0,
                     temp=37.0,
                     noLP=False,
                     circ=False,
                     verb=False):
    """Compute an energy range that computes a given number of structures.

    .. note:: If your RAM is big enough, ``barriers`` can be compiled to read \
      billions of structures, however computations with more than 10.000.000 \
      structures may take a very long time.

    Args:
      seq (str): nucleic acid sequence.
      RNAsubopt (str, optional): path to executable.
      nos (int, optional): number of structures.
      maxe (float, optional): an upper bound of the energy range.
      temp (float, optional): temperature in Celsius.
      noLP (bool): exclude lonely-base-pairs in suboptimal structures. Defaults to False.
      circ (bool): compute density of states. Defaults to False.
      verb (bool): Print verbose information as a comment '#'. Defaults to False.

    Returns:
      [tuple]: (energy-range, number-of-structures)
    """

    if which(RNAsubopt) is None:
        raise ExecError(
            RNAsubopt,
            "RNAsubopt",
            'http://www.tbi.univie.ac.at/RNA')

    num, nump = 0, 0
    e = maxe if nos == 0 else 5.  # Decreasing this value may introduce bugs ...
    ep = e - 1
    while (num < nos + 1):
        if verb:
            print("# Energy: ", "{:.2f}".format(float(e)))

        sptcall = [RNAsubopt, "-D -e {:.2f} -T {:.2f}".format(float(e), temp)]
        if circ:
            sptcall.append("--circ")
        if noLP:
            sptcall.append("--noLP")

        process = sub.Popen([' '.join(sptcall)],
                            stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
        output, err = process.communicate(seq.encode())
        if err:
            # this one might be not so bad ...
            raise Exception(err)

        structures = 0
        for l in output.split(b'\n')[1:-1]:
            [interval, value] = l.split()
            structures += int(value)

            # Make sure that nump fits to the preset ep value
            if nump == 0:
                if ep * 10 == int(interval):
                    nump = structures

            # Break the loop if you reach the number of structures in the
            # output
            if (nos and structures >= nos - (nos * 0.01)
                    ) or float(interval) / 10 >= maxe:
                # print "break", maxe, nos, structures, "E:",
                # float(interval)/10
                e = float(interval) / 10
                num = structures
                break  # return to while

        else:  # if the last loop didn't break, duh!
            num = structures
            if num > nos or num == nump:
                e = float(interval) / 10
                break

            new = e + (e - ep) / (math.log(float(num) / nump)) * \
                (math.log(float(nos) / num))
            ep, nump = e, num
            e = new if new < maxe else maxe
            if abs(e - ep) < 0.1:
                e = ep
                break
            continue
        break  # end while after for loop exit

    return e, num

# ***************** #
# private functions #
# ................. #

def subopt_reaches_minh(fname, minh):
    """ Internal function to report on whether the energy-range of suboptimal
    structures exceeds the ``--minh`` options for ``barriers``. If this is not
    the case, the ``--minh`` option can cause segmentation faults.

    :param fname: Filename of **gzipped and sorted** RNAsubopt result
    :param minh: The ``barriers --minh`` value

    :type fname: string
    :type minh: float

    :return: True or False
    """
    with gzip.open(fname, 'r') as f:
        for i, l in enumerate(f):
            if i == 0:
                continue
            elif i == 1:
                mfe = l.strip().split()[1]
            else:
                sen = l.strip().split()[1]
                if float(sen) - float(mfe) > minh:
                    return True
    return False


def which(program):
    """ Emulates the unix ``which`` command. Snatched from:
      `http://stackoverflow.com/questions/377017/
        test-if-executable-exists-in-python`

      :param program: executable
      :type program: string

      :returns: path-to-executable or None
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def main():
    import ribolands.utils as rnu

    # Quick set test model params """
    name, seq = rnu.parse_vienna_stdin(sys.stdin)
    print(name, seq)

    # param='RNA'
    # dangle='some'
    s_nos = 1e6
    s_maxe = 20
    s_ener = None
    s_patch = []  # ['|', 'pylands_spatch.py', '--theo']

    b_minh = 2.0
    b_maxn = 20

    temp = 37.0
    circ = False
    noLP = True
    verb = True
    force = True
    binrates = True

    if s_ener is None:
        s_ener, s_nos = sys_subopt_range(
            seq, nos=s_nos, maxe=s_maxe, verb=verb)
    sfile = sys_suboptimals(name, seq,
                            ener=s_ener,
                            noLP=noLP,
                            opts=s_patch,
                            verb=verb,
                            force=force)

    [sfile, bfile, efile, rfile, psfile, msfile] = sys_barriers(name, seq, sfile,
                                                        minh=b_minh, maxn=b_maxn, rates=True, binrates=binrates,
                                                        verb=verb, noLP=noLP, force=force)
    tfile, _ = sys_treekin(name, seq, bfile, rfile, binrates=binrates,
                           p0=['2=1'], t0=1e-6, ti=1.02, t8=1e10, verb=verb, force=force)
    pfile = rnu.plot_nxy(name + '.pdf', tfile,
                         ylim=(0, 1), xlim=(1e-2, 1e10), lines=[])

    # BCG = barriersCG(mfile, efile)
    RM = rnu.parse_ratefile(rfile, binary=binrates)
    BT = rnu.parse_barfile(bfile, seq=seq)

    print(RM, BT)
    print(sfile, bfile, efile, rfile, psfile, tfile, pfile)
    return


if __name__ == '__main__':
    main()
