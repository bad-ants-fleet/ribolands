#
#  Coded by: Stefan Badelt <stef@tbi.univie.ac.at>
#  University of Vienna, Department of Theoretical Chemistry
#
#  -*- Content -*-
#  *) systemcalls of RNAsubopt, Kinfold, barriers and treekin
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
from ribolands import (_MIN_VIENNARNA_VERSION, 
                       _MIN_KINFOLD_VERSION,
                       _MIN_BARRIERS_VERSION,
                       _MIN_TREEKIN_VERSION)

# **************************************************************************** #
# Helper functions, Errors, etc.                                               #
# ............................................................................ #

def which(program):
    """Emulates the unix ``which`` command. 
    
    Taken from:
        `http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python`

    Args:
        program (str): executable

    Returns:
        [str]: The path-to-executable or None
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

def check_version(program, rv):
    """Check if the version of an installed program meets minimum requirements.

    Args:
        program (str): executable
        rv (str): requested version number

    Raises:
        VersionError: if version is less than the requested version.

    Returns:
        None: If everything is ok.
    """
 
    def versiontuple(rv):
        return tuple(map(int, (rv.split("."))))

    if which(program) is None:
        if 'treekin' in program:
            raise ExecError(program, "treekin", 
                            'https://www.tbi.univie.ac.at/RNA/Treekin')
        elif 'barries' in program:
            raise ExecError(program, "barriers",
                            'https://www.tbi.univie.ac.at/RNA/Barriers')
        elif 'RNAsubopt' in program:
            raise ExecError(program, "RNAsubopt",
                            'https://www.tbi.univie.ac.at/RNA')
        elif 'Kinfold' in program:
            raise ExecError(program, "Kinfold",
                            'https://www.tbi.univie.ac.at/RNA')
        else:
            raise ExecError(program)

    p, pv = sub.check_output([program, '--version']).split()
    pv = pv.decode()

    if versiontuple(pv) < versiontuple(rv):
        raise VersionError(program, pv, rv)

    return None

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


# **************************************************************************** #
# System call Workflow object                                                  #
# ............................................................................ #

class Workflow(object):
    """ A wrapper class to unify workflows using the syscalls provided.

    Use a Workflow for a specific pipeline, such that input and output files do
    not have to be specified or set manually. If you want to fork a pipeline, e.g.
    to branch off and use different minh for barriers, then use copy() and 
    keep the names of the files up to that point in memory.
    """

    def __init__(self, sequence, model_details, name = 'NoName'):
        # TODOs and other notes ...
        #   - Kinfold does not crash, it reports to STDERR and exits.
        #   - Treekin does not crash with single rate, it reports to STDERR and exits.
        #   - Barriers rate model should allow to adjust the k0 parameter.

        self.md = model_details
        
        # Simple properties, I/O handling:
        self.name = name
        self.sequence = sequence
        self.paramFile = ''
        self.force = False
        self.verbose = False
        self.zip_suboptimals = True
    
        # Prgrams -- need to be set for a mandatory version check.
        self._RNAsubopt = None
        self._barriers = None
        self._treekin = None
        self._Kinfold = None

        # Stuff to calculate
        self.subopt_range = 1 # kcal/mol
        self.subopt_number = None

        # All sorts of files.
        self._outdir = './'
        self.sofile = None # subopt output file
        self.bofile = None # barriers output file
        self.brfile = None # barriers rate file
        self.bbfile = None # barriers binary rate file
        self.bpfile = None # barriers plot file
        self.bmfile = None # barriers mapstruc file
        self.tofile = None # treekin output file
        self.kofile = None # kinfold output file
        self.klfile = None # kinfold log file

        # Some options where only default values are allowed.
        self.moves = 'single-base-pair'
        assert self.md.dangles == 2
        assert self.md.logML == 0
        assert self.md.special_hp == 1
        assert self.md.noGU == 0
        assert self.md.noGUclosure == 0
        assert self.md.gquad == 0
        assert self.md.logML == 0
        assert self.md.circ in [0, 1]
        assert self.md.noLP in [0, 1]
        assert self.md.min_loop_size == 3
        #assert self.md.pf_smooth == 0

    def copy(self):
        other = Workflow(self.sequence, self.md)
        other.name = self.name
        other.paramFile = self.paramFile
        other.force = self.force
        other.verbose = self.verbose
        other.zip_suboptimals = self.zip_suboptimals
    
        # Prgrams -- need to be set for a mandatory version check.
        other._RNAsubopt = self._RNAsubopt
        other._barriers = self._barriers
        other._treekin = self._treekin
        other._Kinfold = self._Kinfold

        # Stuff to calculate
        other.subopt_range = self.subopt_range
        other.subopt_number = self.subopt_number

        # All sorts of files.
        other._outdir = self._outdir
        other.sofile = self.sofile
        other.bofile = self.bofile
        other.brfile = self.brfile
        other.bbfile = self.bbfile
        other.bpfile = self.bpfile
        other.bmfile = self.bmfile
        other.tofile = self.tofile
        other.kofile = self.kofile
        other.klfile = self.klfile
        return other

    def RT(self):
        RT = 0.61632077549999997
        if self.md.temperature != 37.0:
            kelvin = 273.15 + self.md.temperature
            RT = (RT / 310.15) * kelvin
        return RT
    
    @property
    def RNAsubopt(self):
        return self._RNAsubopt

    @RNAsubopt.setter
    def RNAsubopt(self, executable):
        check_version(executable, _MIN_VIENNARNA_VERSION)
        self._RNAsubopt = executable

    @property
    def Kinfold(self):
        return self._Kinfold

    @Kinfold.setter
    def Kinfold(self, executable):
        check_version(executable, _MIN_KINFOLD_VERSION)
        self._Kinfold = executable

    @property
    def barriers(self):
        return self._barriers

    @barriers.setter
    def barriers(self, executable):
        check_version(executable, _MIN_BARRIERS_VERSION)
        self._barriers = executable

    @property
    def treekin(self):
        return self._treekin

    @treekin.setter
    def treekin(self, executable):
        check_version(executable, _MIN_TREEKIN_VERSION)
        self._treekin = executable

    @property
    def outdir(self):
        return self._outdir

    @outdir.setter
    def outdir(self, path):
        if not os.path.exists(path):
            os.mkdir(path)
        self._outdir = path + '/' if path[-1] != '/' else path

    def find_subopt_range(self, nos, maxe, verbose = None):
        if verbose is None:
            verbose = self.verbose

        if self.sequence is None:
            raise Exception('Must set sequence first.')

        if self.RNAsubopt is None:
            raise SystemExit(f'Need to specify executable for RNAsubopt.')

        sener, snum = sys_subopt_range(self.sequence, 
                nos = nos,
                maxe = maxe,
                RNAsubopt = self.RNAsubopt,
                params = self.paramFile,
                temp = self.md.temperature,
                noLP = self.md.noLP,
                circ = self.md.circ,
                verb = verbose)

        self.subopt_range = sener
        self.subopt_number = snum

    def call_RNAsubopt(self, opts = None, force = None, verbose = None):
        if self.RNAsubopt is None:
            raise SystemExit(f'Need to specify executable for RNAsubopt.')

        if force is None:
            force = self.force
        if verbose is None:
            verbose = self.verbose
        
        if True: # more consistent
            sortopt = ['|', 'sort', '-T', '/tmp', '-k3r', '-k2n'] 
        else: # faster
            sortopt = ['-s']

        name = self.outdir + self.name
        sofile = sys_suboptimals(name, self.sequence,
                RNAsubopt = self.RNAsubopt,
                params = self.paramFile,
                ener = self.subopt_range,
                temp = self.md.temperature,
                noLP = self.md.noLP,
                circ = self.md.circ,
                opts = opts,
                sort = sortopt,
                zipped = self.zip_suboptimals,
                force = force,
                verb = verbose)
        
        if verbose and self.sofile:
            print(f'Replacing internal subopt file: {sofile}')

        self.sofile = sofile
        return self.sofile

    def call_barriers(self, 
            rates = True, k0 = 1.0,
            minh = 0.001, maxn = 50, 
            paths = None, bmfile = None, connected = False,
            bsize = False, ssize = False, saddle = False, plot = False, 
            force = None, verbose = None):
        
        if self.barriers is None:
            raise SystemExit(f'Need to specify executable for barriers.')

        if verbose is None:
            verbose = self.verbose
        if force is None:
            force = self.force

        name = self.outdir + self.name
        files = sys_barriers_180(name, self.sofile,
                barriers = self.barriers,
                minh = minh, 
                maxn = maxn,
                temp = self.md.temperature,
                noLP = self.md.noLP,
                circ = self.md.circ,
                moves = self.moves,
                zipped = self.zip_suboptimals,
                rates = rates,
                k0 = k0,
                paths = paths,
                bsize = bsize,
                ssize = ssize,
                saddle = saddle,
                bmfile = bmfile,
                connected = connected,
                plot = plot,
                force = force,
                verbose = verbose)

        if verbose and self.bofile:
            print(f'Replacing internal barriers output file: {files[0]}')
        if verbose and self.brfile:
            print(f'Replacing internal barriers rates file: {files[2]}')
        if verbose and self.bbfile:
            print(f'Replacing internal barriers binary rates file: {files[3]}')
        if verbose and self.bpfile:
            print(f'Replacing internal barriers plot file: {files[4]}')
        if verbose and self.bmfile:
            print(f'Replacing internal barriers mapstruc file: {files[5]}')

        self.bofile = files[0]
        # efile 
        self.brfile = files[2]
        self.bbfile = files[3]
        self.bpfile = files[4]
        self.bmfile = files[5]

        return files

    def call_treekin(self, 
            p0 = ['1=0.5', '2=0.5'], 
            t0 = 0, ti = 1.2, t8 = 1e6,
            binrates = True, bofile = None, # bofile = False to avid using internal one.
            useplusI = False, exponent = False, mpack_method = None,
            quiet = True, force = None, verbose = None):

        if self.treekin is None:
            raise SystemExit(f'Need to specify executable for treekin.')

        if force is None:
            force = self.force
        if verbose is None:
            verbose = self.verbose

        ratefile = self.bbfile if binrates else self.brfile
        if bofile is False:
            bofile is None
        elif not bofile:
            bofile = self.bofile

        name = self.outdir + self.name
        tofile, tefile = sys_treekin_051(name, ratefile,
                treekin = self.treekin,
                bofile = bofile,
                p0 = p0,
                t0 = t0,
                ti = ti,
                t8 = t8,
                binrates = binrates,
                useplusI = useplusI,
                exponent = exponent,
                mpack_method = mpack_method,
                quiet = quiet,
                force = force,
                verbose = verbose)

        if verbose and self.tofile:
            print(f'Replacing internal treekin output file: {tofile}')

        self.tofile = tofile
        return tofile, tefile

    def call_Kinfold(self, start, stop, fpt = True, rect = False,
            time = 5000, num = 1, 
            ratemodel = 'Metropolis',
            lmin = False, silent = True, erange = 100,
            name = None,
            force = None, verbose = None):

        if self.Kinfold is None:
            raise SystemExit(f'Need to specify executable for Kinfold.')

        if verbose is None:
            verbose = self.verbose

        if name is None: 
            name = self.outdir + self.name
        klfile, kefile, kofile = sys_kinfold(name, self.sequence, 
                kinfold = self.Kinfold,
                start = start,
                stop = stop,
                fpt  = fpt,
                rect = rect,
                time = time,
                num = num,
                ratemodel = ratemodel,
                moves = self.moves,
                noLP = self.md.noLP,
                logML = self.md.logML,
                dangle = self.md.dangles,
                temp = self.md.temperature,
                params = self.paramFile,
                erange = erange,
                lmin = lmin,
                silent = silent,
                force = force,
                verbose = verbose)

        if verbose and self.klfile:
            print(f'Replacing internal kinfold log file: {klfile}')
        if verbose and self.kofile:
            print(f'Replacing internal kinfold output file: {kofile}')

        self.klfile = klfile
        self.kofile = kofile
        return klfile, kefile, kofile


# **************************************************************************** #
# System call wrapper functions                                                #
# ............................................................................ #

#Last update: version 0.5.1
def sys_treekin_051(basename, ratefile,
                treekin = 'treekin',
                bofile = None,
                p0 = ['1=0.5', '2=0.5'],
                t0 = 0,
                ti = 1.2,
                t8 = 1e6,
                binrates = True,
                useplusI = False,
                exponent = False,
                mpack_method = None,
                quiet = True,
                force = False,
                verbose = False):
    """Perform a system-call of the program ``treekin``.

    Prints the results into files and returns the respective filenames. This
    wrapper will produce two files, one from ``STDIN`` and ``STDERR`` output.

    Args:
        basename (str): Basename of output files.
        ratefile (str): Input rate matrix. 
            Fileformat as produced by ``barriers``.
        treekin (str, optional): Path to executable. Defaults to treekin.
        bofile (str, optional): Optional barriers file. Defaults to None.
        p0 ([str,...], optional): A list to specify an initial occupancy vector. 
            Defaults to: ['1=0.5', '2=0.5']
        t0 (float, optional): Start time of simulation.
        ti (float, optional): Time increment of the treekin solver.
        t8 (float, optional): Stop time for the solver.
        binrates (bool, optional): Ratefile is in binary format. 
            Defaults to True.
        useplusI (bool, optional): Use treekin method --useplusI. 
            Defaults to False.
        exponent (bool, optional): Use treekin method --exponent. 
            Defaults to False.
        mpack_method (str, optional): Set treekin parameter --mpack-method. 
            Defaults to None.
        force (bool, optional): Overwrite existing files. Defaults to False.
        verbose (bool, optional): Print verbose information. Defaults to False.

    Raises:
        ExecError: Program does not exist.
        SuboprocessError: Program terminated with exit code: ...

    Returns:
        [str, str]: The name of the file containing STDOUT and STDERR of
            ``treekin`` call.  
    """

    if which(treekin) is None:
        raise ExecError(treekin,
                        "treekin", 'http://www.tbi.univie.ac.at/RNA/Treekin')

    reg_flt = re.compile(b'[-+]?[0-9]*.?[0-9]+([eE][-+]?[0-9]+)?.')
    # http://www.regular-expressions.info/floatingpoint.html

    tofile = basename + '.tkn'
    tefile = basename + '.err'

    if not force and os.path.exists(tofile):
        if verbose:
            print("# {:s} <= Files exist".format(tofile))
        return tofile, tefile

    # Unfortunately, running treekin with a single state leads to an error that
    # is printed to STDOUT instead of STDERR. The program, exits with success.
    # That's why we catch it first:
    if binrates:
        with open(ratefile, 'rb') as rf:
            i, = unpack('i', rf.read(calcsize('i')))  # unsigned long, little-endian
        i -= 1
    else:
        with open(ratefile) as rf:
            for i, _ in enumerate(rf):
                pass
    if i == 0:
        raise SubprocessError(None, 'No transition rates found.')

    syscall = [treekin, '--method', 'I']
    if quiet:
        syscall.extend(['--quiet'])
    if useplusI:
        syscall.extend(['--useplusI'])
    if mpack_method:
        syscall.extend(['--mpack-method=' + mpack_method])
    syscall.extend(['--tinc', str(ti)])
    syscall.extend(['--t0', str(t0)])
    syscall.extend(('--t8', str(t8)))
    for p in p0:
        syscall.extend(('--p0', p))
    if bofile:
        syscall.extend(('--bar', bofile))
    if binrates:
        syscall.extend(['--bin'])
    if exponent:
        syscall.extend(['--exponent'])

    call = "cat {} | {} 2> {} > {}".format(ratefile, ' '.join(syscall), tefile, tofile)

    if verbose:
        print(f'# {call}')

    # Do the simulation (catch treekin errors)
    with open(ratefile, 'r') as rts, \
            open(tofile, 'w') as tkn, \
            open(tefile, 'w') as err:
        proc = sub.Popen(syscall, stdin = rts, stdout = tkn, stderr = err)
        proc.communicate(None)
        if proc.returncode:
            raise SubprocessError(proc.returncode, call)

    return tofile, tefile

#Last update: version 1.8.0
def sys_barriers_180(basename, sofile,
                 barriers = 'barriers',
                 minh = 0.001,
                 maxn = 50,
                 temp = 37.0,
                 noLP = False,
                 moves = 'single-base-pair',
                 zipped = True,
                 rates = True,
                 k0 = 1.,
                 paths = None,
                 bsize = False,
                 ssize = False,
                 circ = False,
                 saddle = False,
                 bmfile = None,
                 plot = False,
                 connected = False,
                 force = False,
                 verbose = False):
    """Perform a system-call of the program ``barriers``.

    The print the results into a file and return the respective filename. This
    wrapper will return the output files, including ``STDIN`` and ``STDERR``.

    Args:
      name (str): Name of the sequence used to name the output file.
      sofile (str): Input filename for ``barriers`` as produced by ``RNAsubopt``.
      temp (float, optional): Specify temperature in Celsius.
      force (bool, optional): Overwrite existing files with the same name.
      verbose (bool, optional): Print the current system-call to ``stderr``.

    Raises:
      ExecError: Program does not exist.
      SuboprocessError: Program terminated with exit code: ...

    Returns:
      [bofile, befile, brfile, bbfile, bpfile, msfile]: 
        A list of produced files containing ``barriers`` results.
    """

    if zipped and which('zcat') is None:
        print('Using gzipped subopt files requires the commandline tool zcat.')
        raise ExecError('zcat', "zcat")

    if which(barriers) is None:
        raise ExecError(barriers, "barriers",
                        'http://www.tbi.univie.ac.at/RNA/Barriers')

    if not sofile or not os.path.exists(sofile):
        raise Exception('Cannot find input file:', sofile)

    bofile = basename + '.bar'
    befile = basename + '.err'
    brfile = basename + '_rates.txt'
    bbfile = basename + '_rates.bin'
    bpfile = basename + '_tree.ps' if plot else None
    msfile = basename + '.ms' if bmfile else None

    if not force and os.path.exists(bofile) and os.path.exists(brfile) and os.path.exists(bbfile) and \
      (not paths or all(map(os.path.exists, ['{:s}_path.{:03d}.{:03d}.txt'.format(basename, int(x), int(y)) for x, y in map(lambda x: x.split('='), paths)]))) and \
      (not plot or os.path.exists(bpfile)) and \
      (not bmfile or os.path.exists(msfile)):
        if verbose:
            print("#", bofile, brfile, " <= Files exist")
        return [bofile, befile, brfile, bbfile, bpfile, msfile]

    barcall = [barriers]

    if not plot:
        barcall.extend(['-q'])

    if connected:
        barcall.extend(['-c'])

    if paths:
        for p in paths:
            barcall.extend(['--path', p])

    if noLP:
        barcall.extend(['-G', 'RNA-noLP'])
    else:
        barcall.extend(['-G', 'RNA'])

    if moves == 'single-base-pair':
        pass
    elif moves == 'shift':
        barcall.extend(['--moves=Shift'])
    else:
        raise ValueError(f"Invalid move-set for barriers: {moves}")

    # buggy barriers
    if subopt_reaches_minh(sofile, minh, zipped):
        barcall.extend(["--minh", str(minh)])
        barcall.extend(["--max", str(int(maxn))])
    barcall.extend(["-T", str(temp)])

    if rates:
        barcall.extend(['--rates'])
        barcall.extend(['--rates-text-file', brfile])
        barcall.extend(['--rates-binary-file', bbfile])
    if bsize:
        barcall.extend(['--bsize'])
    if ssize:
        barcall.extend(['--ssize'])
    if saddle:
        barcall.extend(['--saddle'])

    if bmfile:
        barcall.extend(["--mapstruc", bmfile])
        barcall.extend(["--mapstruc-output", msfile])

    call = "{} 2> {} > {}".format(' '.join(barcall), befile, bofile)
    if verbose:
        print(f'# {call}')

    # TODO: This version is prettier than the old one, but it might be slower
    # than just writing the string and "shell = True".
    if zipped:
        inp = sub.Popen(["zcat", sofile], stdout = sub.PIPE)
        with open(bofile, 'w') as bh, open(befile, 'w') as eh:
            proc = sub.Popen(barcall, stdin = sub.PIPE, stdout = bh, stderr = eh)
            proc.communicate(inp.communicate()[0])
            if proc.returncode:
                raise SubprocessError(proc.returncode, call)
    else:
        with open(sofile, 'r') as sh, open(bofile, 'w') as bh, open(befile, 'w') as eh:
            proc = sub.Popen(barcall, stdin = sh, stdout = bh, stderr = eh)
            proc.communicate()
            if proc.returncode:
                raise SubprocessError(proc.returncode, call)

    if rates and k0 != 1.: # So this should actually be supported by barriers, but it's not.
        with open(bbfile, 'rb') as rf, \
                open(brfile + '.tmp', 'w') as nr, \
                open(bbfile + '.tmp', 'wb') as nb:
            dim, = unpack('i', rf.read(calcsize('i')))
            nb.write(pack("i", dim))

            rm = []
            for e in range(dim):
                col = []
                for e in range(dim):
                    r, = unpack('d', rf.read(8))
                    rate = r * k0
                    nb.write(pack("d", rate))
                    col.append(rate)
                rm.append(col)

            for line in zip(*rm):
                newline = "".join(map("{:10.4g}".format, line))
                nr.write(newline + "\n")
 
        os.rename(brfile + '.tmp', brfile)
        os.rename(bbfile + '.tmp', bbfile)

    if plot:
        os.rename('tree.ps', bpfile)

    if paths:
        for p in paths:
            x, y = p.split('=')
            pfname = 'path.{:03d}.{:03d}.txt'.format(int(x), int(y))
            os.rename(pfname, basename + '_' + pfname)

    return [bofile, befile, brfile, bbfile, bpfile, msfile]


def sys_suboptimals(name, seq,
                    RNAsubopt = 'RNAsubopt',
                    params = None,
                    ener = None,
                    temp = 37.0,
                    noLP = False,
                    circ = False,
                    opts = None,
                    sort = ['|', 'sort', '-T', '/tmp', '-k3r', '-k2n'],
                    zipped = True,
                    force = False,
                    verb = False):
    """Perform a system-call of the program ``RNAsubopt``.

    The print the results into a file and return the filename.

    Args:
        name (str): Name of the sequence used to name the output file.
        seq (str): Nucleic acid sequence.
        RNAsubopt (str, optional): Path to ``RNAsubopt'' executable.
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

    Returns:
        [str]: Filename of the file containing ``RNAsubopt`` results.
    """

    if which(RNAsubopt) is None:
        raise ExecError(RNAsubopt, "RNAsubopt", 'http://www.tbi.univie.ac.at/RNA')

    if ener is None:
        ener, nos = sys_subopt_range(seq, verb=verb,
                                     RNAsubopt=RNAsubopt, noLP=noLP, circ=circ, temp=temp)
        if verb:
            print("# Energy-Update: {:.2f} kcal/mol to compute {} sequences".format(
                ener, nos))

    sofile = name + '.spt.gz' if zipped else name + '.spt'

    if not force and os.path.exists(sofile):
        if verb:
            print("#", sofile, "<= File exists")
        return sofile

    sptcall = [RNAsubopt, "-e {:.2f} -T {:.2f}".format(ener, temp)]
    if params:
        sptcall.extend(["--paramFile", params])
    if noLP:
        sptcall.append("--noLP")
    if circ:
        sptcall.append("--circ")
    if opts:
        sptcall.extend(opts)
    sptcall.extend(sort)

    if zipped:
        sptcall.extend(('|', 'gzip', '--best'))

    if verb:
        print("#", "echo \"{}\" | {} > {}".format(seq, ' '.join(sptcall), sofile))

    with open(sofile, 'w') as shandle:
        proc = sub.Popen([' '.join(sptcall)],
                         stdin=sub.PIPE, stdout=shandle, shell=True)
        proc.communicate(seq.encode())
        if proc.returncode:
            call = "echo \"{}\" | {} > {}".format(seq, ' '.join(sptcall), sofile)
            raise SubprocessError(proc.returncode, call)

    return sofile


def sys_subopt_range(seq,
                     RNAsubopt='RNAsubopt',
                     params = None,
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
        raise ExecError(RNAsubopt, "RNAsubopt", 'http://www.tbi.univie.ac.at/RNA')

    num, nump = 0, 0
    e = maxe if nos == 0 else 5.  # Decreasing this value may introduce bugs ...
    ep = e - 1
    while (num < nos + 1):
        if verb:
            print("# Energy: ", "{:.2f}".format(float(e)))

        sptcall = [RNAsubopt, "-D -e {:.2f} -T {:.2f}".format(float(e), temp)]
        if params:
            sptcall.extend(["--paramFile", params])
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

            # Break the loop if you reach the number of structures in the output
            if (nos and structures >= nos - (nos * 0.01)) or float(interval) / 10 >= maxe:
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

def sys_kinfold(name, seq, 
        kinfold = 'Kinfold',
        start = None,
        stop = None,
        fpt  = False,
        rect = False,
        time = 5000,
        num = 1,
        ratemodel = 'Metropolis',
        moves = 'single-base-pair',
        noLP = False,
        logML = False,
        dangle = 2,
        temp = 37,
        params = None,
        # Output
        erange = 20, # Kinfold Default.
        lmin = False,
        silent = False,
        force = False,
        verbose = False):
    """Perform a system-call of the program ``Kinfold``.

    The print the results into a file and return the respective filename. This
    wrapper will return the output files, including ``STDIN`` and ``STDERR``.

    Args:
      ...

    Returns:
      ...
    """

    name += '_kinfold'

    klfile = name + '.log'
    kefile = name + '.err'
    kofile = name + '.out'

    if os.path.exists(klfile) and os.path.exists(kofile):
        if force is None:
            if verbose:
                print("#", klfile, kofile, "<= Files exist, appending output!")
        elif force is False:
            return [klfile, kefile, kofile]
        else:
            if verbose:
                print("#", klfile, kofile, "<= Removing old files!")
            if os.path.exists(klfile):
                os.remove(klfile)
            if os.path.exists(kofile):
                os.remove(kofile)

    kinput = seq + '\n'

    syscall = [kinfold]
    syscall.extend(['--num', str(int(num))])
    syscall.extend(['--time', str(time)])
    syscall.extend(['--log', name])
    syscall.extend(['--cut', str(erange)])
    #syscall.extend(['--seed', "62159=58010=26254"])
    if params:
        syscall.extend(["--Par", params])
    if dangle != 2:
        syscall.extend(['--dangle', str(dangle)])
    if temp != 37:
        syscall.extend(['-T', str(temp)])
    if ratemodel == 'Kawasaki':
        pass
    elif ratemodel == 'Metropolis':
        syscall.extend(['--met'])
    else:
        raise NotImplementedError('unknown rate model')

    if lmin:
        syscall.extend(['--lmin'])

    if not fpt: # NOTE: fpt switches first passage time off (!!!!)
        syscall.extend(['--fpt'])

    if rect: # NOTE: rect is only supported through a development branch right now ...
        syscall.extend(['--rect'])

    if not logML: # NOTE: logML switches logarithmic multiloop evaluation off (!!!!)
        syscall.extend(['--logML'])

    if silent:
        syscall.extend(['--silent'])

    if noLP:
        syscall.extend(['--noLP'])

    if moves == 'single-base-pair':
        syscall.extend(['--noShift'])
    elif moves == 'shift':
        pass

    if start:
        assert isinstance(start, str)
        syscall.extend(['--start'])
        kinput += start + '\n'

    if stop:
        syscall.extend(['--stop'])
        if isinstance(stop, str):
            kinput += stop + '\n'
        else:
            assert isinstance(stop, list)
            kinput += '\n'.join(stop)

    if verbose:
        print('#', ' '.join(syscall), '2>', kefile, '>', kofile)

    with open(kofile, 'w') as khandle, open(kefile, 'w') as ehandle:
        proc = sub.Popen(' '.join(syscall), 
                stdin=sub.PIPE, stdout=khandle, stderr=ehandle, shell=True)
        proc.communicate(kinput.encode())
        if proc.returncode:
            call = "{} 2> {} > {}".format(
                ' '.join(syscall), kefile, kofile)
            print(proc.returncode)
            print(call)

    return [klfile, kefile, kofile]

# ***************** #
# private functions #
# ................. #

def subopt_reaches_minh(fname, minh, zipped = True):
    """ Internal function to report on whether the energy-range of suboptimal
    structures exceeds the ``--minh`` options for ``barriers``. If this is not
    the case, the ``--minh`` option can cause segmentation faults.

    :param fname: Filename of **gzipped and sorted** RNAsubopt result
    :param minh: The ``barriers --minh`` value

    :type fname: string
    :type minh: float

    :return: True or False
    """
    def searchit(f):
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

    if zipped:
        with gzip.open(fname, 'r') as f:
            return searchit(f)
    else:
        with open(fname, 'r') as f:
            return searchit(f)

    raise AssertionError('Function must exit above.')

# ===========
# DEPRECATED
# ===========
def sys_treekin(name, seq, bofile, brfile,
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
    print('WARNING: Using deprecated function: sys_treekin, use sys_treekin_051 instead.')

    return sys_treekin_051(name, brfile,
                treekin = treekin,
                bofile = bofile,
                p0 = p0,
                t0 = t0,
                ti = ti,
                t8 = t8,
                binrates = binrates,
                useplusI = useplusI,
                exponent = exponent,
                mpack_method = mpack_method,
                force = force,
                verbose = verb)


def sys_barriers(name, seq, sofile,
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
                 ssize=False,
                 circ=False,
                 saddle=False,
                 bmfile='',
                 force=False,
                 verb=False):

    print('WARNING: Using deprecated function: sys_barriers, use sys_barriers_180 instead.')

    [bofile, befile, brfile, bbfile, bpfile, msfile] = sys_barriers_180(name, sofile,
                 barriers = barriers,
                 minh = minh, 
                 maxn = maxn,
                 temp = temp,
                 noLP = noLP,
                 moves = moves,
                 zipped = gzip,
                 rates = rates,
                 k0 = k0,
                 bsize = bsize,
                 ssize = ssize,
                 circ = circ,
                 saddle = saddle,
                 bmfile = bmfile,
                 force = force,
                 verbose = verb)
    if binrates:
        return [None, bofile, befile, bbfile, None, msfile]
    else:
        return [None, bofile, befile, brfile, None, msfile]
 
def main():
    import RNA 
    from ribolands.utils import parse_vienna_stdin, parse_ratefile, plot_nxy, plot_nxy_linlog
    from ribolands.parser import parse_barriers

    #if paramfile:
    #    RNA.read_parameter_file(args.paramFile)

    # Set model details.
    vrna_md = RNA.md()
    vrna_md.noLP = 0
    vrna_md.temperature = 25
    vrna_md.dangles = 2
    vrna_md.logML = 0
    vrna_md.special_hp = 1
    vrna_md.noGU = 0
    
    # Quick set test model params """
    name, seq = parse_vienna_stdin(sys.stdin)

    Pipe = Workflow(seq, vrna_md, name = name)
    Pipe.force = True
    Pipe.verbose = True
    Pipe.outdir = 'mypipe'
    Pipe.RNAsubopt = 'RNAsubopt'
    Pipe.barriers = 'barriers'
    Pipe.treekin = 'treekin'
    Pipe.Kinfold = 'Kinfold'
    Pipe.zip_suboptimals = False

    Pipe.find_subopt_range(nos = 10000, maxe = 20)

    sofile = Pipe.call_RNAsubopt()
    bofile, befile, brfile, bbfile, *_ = Pipe.call_barriers(minh = 3, maxn = 20, plot = False, 
        connected = True, force = True, verbose = True)

    tofile, _ = Pipe.call_treekin(p0 = ['2=1'], useplusI = True, verbose = True)


    lmins = parse_barriers(bofile, return_tuple = True) 
    assert lmins[0] == Pipe.sequence

    for lm in lmins:
        print(lm)

    stop = []
    for lm in lmins[1:-1]:
        stop.append(lm.structure)

    klfile, *_ = Pipe.call_Kinfold(start = lmins[-1].structure, stop = stop, time = 1e5, num = 5)

    pfile = plot_nxy(Pipe.name + '.pdf', tofile, ylim = None, xlim = (1, 1e5), lines = [])

    pfile = plot_nxy_linlog(Pipe.name + '_linlog.pdf', tofile, xlim = (0, 1000, 1e5), figdivide = 1.5)

    RM = parse_ratefile(bbfile, binary = True)

    for l in RM:
        print(l)

if __name__ == '__main__':
    main()
