#!/usr/bin/env python

"""
  BarMap 2 -- cotranscriptional folding with *barriers* and *treekin*

  :requires: python-v.2.7, RNAsubopt, barriers-v1.6, treekin-v0.4

  use two whitespace characters as tab when editing this file!
  vim-config = set: ts=2 et sw=2 sts=2
"""

import os
import re
import sys
import argparse
import numpy as np
import subprocess as s  # barmap_treekin()
import collections as c
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import ribolands as ril
from ribolands.syswraps import which, ExecError, SubprocessError


class LostPopulationError(Exception):
    """Raise Error: Lost significant population."""

    def __init__(self, var):
        self.message = "{}".format(var)
        super(LostPopulationError, self).__init__(self.message)


def get_plot_data(tfiles, plist, args):
    """

    """
    verb = args.verbose
    start = args.start
    stop = args.stop
    t8 = args.t8
    tX = args.tX
    plim = 0.02  # former argument plot-cutoff, seems unnecessary
    # NOTE: setting plim to a very small number leads to vertical lines in the plot. Why?
    # Because sometimes low populated minima merge into a high populated minimum, and then
    # this will be visible.

    tot_time = 0.
    # List of trajectories
    tlist = [[] for i in range(len(plist) + 1)]

    reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')
    # http://www.regular-expressions.info/floatingpoint.html

    for el, l in enumerate(range(start, stop + 1)):
        tfile = tfiles[el]

        # Get raw data => nxy
        #  t 1 2 3 4
        # [x y y y y]
        # [x y y y y]
        nxy = []
        with open(tfile) as tkn:
            for line in tkn:
                if re.match('#', line):
                    continue
                elif not reg_flt.match(line):
                    raise Exception(tfile, "neither comment nor time course!")
                nxy.append(map(float, line.strip().split()))
        if nxy == []:
            nxy.append([t8, 1.0])

        # Make list of trajectories
        # t [x x x x x]
        # 1 [y y y y y]
        # 2 [y y y y y]
        # 3 [y y y y y]
        # 4 [y y y y y]
        #
        # and use the pathlist with lmin information
        # [0 0 0 3 3 2 2 2 0]
        # [1 1 1 5 5 6 8 9 0]
        # [3 3 3 2 2 1 1 2 1]
        # [2 2 2 2 2 1 1 2 1]
        # [0 9 9 9 9 9 1 2 1]
        # ...
        #
        # to replace every lmin at the current transcription step with the
        # trajectory (or None) => avoid duplicates with seen=set()!

        traject = map(list, zip(*nxy))
        timeline = map(lambda x: x + tot_time, traject[0])
        tot_time = timeline[-1]

        tlist[0].extend([timeline])

        # replace every position in tlist with the time-course of the treekin
        # output (or None). Use the mapping provided in the pathlist for the the
        # current transcription step l

        # Remove seen from this part to get all trajectries from start to end.
        seen = set()
        pmapping = zip(*plist)[l - start]
        for e, idx in enumerate(pmapping, 1):
            # print start, l, pmapping, idx # traject, idx
            # Make sure that we are above the limit of detection
            if (idx == 0) or (idx in seen):
                line = [None]  # for i in range(len(timeline))]
            else:
                seen.add(idx)
                if max(traject[idx]) > plim:
                    line = traject[idx]
                elif el > 1 and tlist[e][-1][-1] is not None:
                    line = traject[idx]
                else:
                    line = [None]  # for i in range(len(timeline))]
            tlist[e].extend([line])

    return tlist


def plot_xmgrace(all_in, plist, args):
    head = """
@with line
@line on
@line loctype world
@line g0
@line linewidth .1
@line linestyle 1
@line color 7
@line arrow 0
@line arrow type 0
@line arrow length 1.000000
@line arrow layout 1.000000, 1.000000
@line def

@map color 0 to (255, 255, 255), "white"
@map color 1 to (0, 0, 0), "black"
@map color 2 to (255, 0, 0), "red"
@map color 3 to (0, 255, 0), "green"
@map color 4 to (0, 0, 255), "blue"
@map color 5 to (255, 255, 0), "yellow"
@map color 6 to (188, 143, 143), "brown"
@map color 7 to (220, 220, 220), "grey"
@map color 8 to (148, 0, 211), "violet"
@map color 9 to (0, 255, 255), "cyan"
@map color 10 to (255, 0, 255), "magenta"
@map color 11 to (255, 165, 0), "orange"
@map color 12 to (114, 33, 188), "indigo"
@map color 13 to (103, 7, 72), "maroon"
@map color 14 to (64, 224, 208), "turquoise"
@map color 15 to (0, 139, 0), "green4"

"""
    best = range(1, 16)
    gfile = args.name + '.gr'

    c = 0
    with open(gfile, 'w') as gfh:
        gfh.write(head)
        for e, course in enumerate(all_in):
            if e == 0:
                continue
            finalmin = plist[e - 1][-1]
            color = best[finalmin - 1] if len(best) > finalmin else best[-1]

            flag = 0
            for x, y in zip(all_in[0], all_in[e]):
                if y[0] is None:
                    continue
                if flag == 0:
                    gfh.write("@  s{} line color {}\n".format(c, color))
                    if color == best[-1]:
                        gfh.write("@  s{} line linewidth 1\n".format(c))
                    else:
                        gfh.write("@  s{} line linewidth 2\n".format(c))
                for i in range(len(x)):
                    gfh.write("{:f} {:f}\n".format(x[i], y[i]))
                flag = 1

            if flag:
                gfh.write("&\n")
                c += 1

    return gfile


def plot_matplotlib(name, seq, tlist, plist, args):
    """ Description """
    start = args.start
    stop = args.stop

    t8 = args.t8
    lin_time = (stop - start) * float(t8)
    tX = lin_time + args.tX if (lin_time +
                                args.tX) >= lin_time * 10 else lin_time * 10

    name = args.name
    title = args.name

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax.set_ylim([0, 1.01])
    ax.set_xscale('linear')

    # Make the second part of the plot logarithmic
    ax.set_xlim((0, lin_time))
    divider = make_axes_locatable(ax)
    axLog = divider.append_axes("right", size=2.5, pad=0, sharey=ax)
    axLog.set_xscale('log')
    axLog.set_xlim((lin_time + 0.00001, tX))
    axLog.set_ylim([0, 1.01])
    axLog.yaxis.set_visible(False)

    best = ['black', 'red', 'green', 'blue', 'yellow', 'brown', 'grey', 'violet',
            'magenta', 'orange', 'indigo', 'maroon', 'cyan']
    seen = set()
    for e, traject in enumerate(tlist):
        if e == 0:
            for i in range(len(traject)):
                ax.axvline(x=traject[i][-1], linewidth=0.1,
                           color='black', linestyle='--')
            continue

        fulltime = []
        fulltrajectory = []
        for i in range(len(traject)):
            if traject[i][0] is None:
                continue
            fulltime += tlist[0][i]
            fulltrajectory += traject[i]

        finalmin = plist[e - 1][-1]

        color = best[finalmin - 1] if len(best) > finalmin else best[-1]
        lin, = ax.plot(fulltime, fulltrajectory, '-', color=color)
        log, = axLog.plot(fulltime, fulltrajectory, color=color)

        if finalmin not in seen:
            seen.add(finalmin)
            # NOTE: Adjust here for legend
            if traject[-1] and max(traject[-1]) > 0.1:
                # if fulltrajectory and max(fulltrajectory) > 0.1:
                log.set_label("lmin {:d}".format(finalmin))

    fig.set_size_inches(7, 3)
    fig.text(0.5, 0.95, title, ha='center', va='center')

    plt.legend()

    ax.set_ylabel('occupancy [mol/l]', fontsize=11)
    ax.set_xlabel('time [seconds]', ha='center', va='center', fontsize=11)
    ax.xaxis.set_label_coords(.9, -0.15)

    # plt.show()
    pfile = name + '.pdf'
    plt.savefig(pfile, bbox_inches='tight')

    return pfile


def barmap_treekin(bname, seq, bfiles, plist, args):
    """ write treekin files from bfiles and plist information """
    tmpdir = args.tmpdir
    start = args.start
    stop = args.stop
    verb = args.verbose
    tX = args.tX
    t8 = args.t8
    cutoff = args.occupancy_cutoff

    p0 = args.p0
    tt = 0
    tfiles = []
    reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')
    # http://www.regular-expressions.info/floatingpoint.html

    for e, l in enumerate(range(start, stop + 1)):
        t8 = tX if l == stop else t8

        cseq = seq[0:l]
        cname = "{}-t8_{}-len_{}".format(bname, t8, l)
        [bfile, efile, rfile, psfile] = bfiles[e]

        with open(bfile) as bf:
            for i, _ in enumerate(bf):
                pass
        if i >= 2:
            try:
                ctfile, _ = ril.sys_treekin(cname, cseq, bfile, rfile,
                                            treekin=args.treekin,
                                            exponent=False,
                                            useplusI=False,
                                            binrates=True,
                                            p0=p0,
                                            t0=args.t0,
                                            ti=args.ti,
                                            t8=t8,
                                            verb=verb,
                                            force=args.force)
            except SubprocessError:
                print "# repeating treekin calculations with --exponent"
                ctfile, _ = ril.sys_treekin(cname, cseq, bfile, rfile,
                                            treekin=args.treekin,
                                            exponent=True,
                                            useplusI=False,
                                            binrates=True,
                                            p0=p0,
                                            t0=args.t0,
                                            ti=args.ti,
                                            t8=t8,
                                            verb=verb,
                                            force=True)

            lastlines = s.check_output(
                ['tail', '-2', ctfile]).strip().split("\n")
            if reg_flt.match(lastlines[0]):
                tt += float(lastlines[0].split()[0])
                if l < stop:
                    curlmin = np.array(plist)[:, e]
                    newlmin = np.array(plist)[:, e + 1]
                    p0 = set_p0(bfile, l, lastlines, curlmin, newlmin,
                                cutoff, verb)
                else:
                    for (i, pop) in enumerate(lastlines[0].split()):
                        if i != 0 and float(pop) > cutoff:
                            ss, en = get_structure(bfile, i, energy=True)
                            print "{:3d} {:3d} {:f} {:s} {:6.2f}".format(l, i,
                                                                         float(pop), ss, float(en))
            else:
                raise SubprocessError(
                    'found crashed treekin trajectory.', ctfile)
        else:
            # Create an empty file
            ctfile = cname + '.tkn'
            open(ctfile, 'a').close()
            tt += t8
            print "{:3d} {:3d} {:f} {:s} {:6.2f}".format(l, 1, 1.0, get_structure(bfile, 1), 0.00)
        tfiles.append(ctfile)
    return tfiles


def barmap_mapping(_bname, seq, args):
    """ Parse mapping info into pathlist """
    tmpdir = args.tmpdir
    start = args.start
    stop = args.stop
    verb = args.verbose
    force = args.force

    plist = []
    mapinfo = "{}-map_{:d}_to_{:d}.map".format(_bname, start, stop)

    if os.path.exists(mapinfo) and not force:
        if verb:
            print "# {} <= File exists".format(mapinfo)
        with open(mapinfo, 'r') as m:
            m.readline().strip()
            for line in m:
                path = [int(x) for x in line.strip().split(" => ")]
                plist.append(path)
    else:
        mlist = []
        for l in range(start, stop + 1):
            cseq = seq[0:l]
            cname = "{}-len_{:02d}".format(_bname, l)
            pname = "{}-len_{:02d}".format(_bname, l - 1)

            if os.path.exists(pname + '.bar'):
                if verb:
                    print "# Get mapping info {:d} -> {:d}".format(l - 1, l)
                mlist.append(get_mapping_dict(pname + '.bar', cname + '.err'))

        plist = pathlist(mlist)
        with open(mapinfo, 'w') as m:
            m.write("# Mapping information for barrier trees {:d} to {:d}\n".format(
                start, stop))
            for path in plist:
                m.write(" => ".join(map("{:4d}".format, path)) + "\n")
    return plist


def barmap_barriers(_bname, seq, sfiles, args):
    """ Print barriers files and mapping info """

    bfiles = []
    prog = ril.ProgressBar(args.stop - args.start + 1)
    for e, l in enumerate(range(args.start, args.stop + 1)):
        cseq = seq[0:l]
        cname = "{}-len_{:02d}".format(_bname, l)
        pname = "{}-len_{:02d}".format(_bname, l - 1)
        sfile = sfiles[e]

        if os.path.exists(pname + '.bar'):
            with open(pname + '.bar', 'r') as oldbar, \
                    open(pname + '.map', 'w') as mapfile:
                for e, line in enumerate(oldbar):
                    if e == 0:
                        continue
                    cols = line.strip().split()
                    mapfile.write(cols[1] + '.' + "\n")
            mfile = pname + '.map'
        else:
            mfile = ''

        # Make sure the first round for mapping is always recomputed
        # force = True if e == 1 else args.force

        [sfile, bfile, efile, rfile, psfile] = ril.sys_barriers(cname, cseq, sfile,
                                                                barriers=args.barriers,
                                                                minh=args.b_minh,
                                                                maxn=args.b_maxn,
                                                                k0=args.k0,
                                                                temp=args.temperature,
                                                                noLP=args.noLP,
                                                                moves='single-base-pair',
                                                                gzip=True,
                                                                rates=True,
                                                                binrates=True,
                                                                bsize=False,
                                                                saddle=False,
                                                                mfile=mfile,
                                                                force=args.force,
                                                                verb=args.verbose)
        bfiles.append([bfile, efile, rfile, psfile])
        prog.inc()
    return bfiles


def barmap_subopts(_sname, seq, args):
    """ Print all suboptimal structures, return list of files """
    sfiles = []
    prog = ril.ProgressBar(args.stop - args.start + 1)
    for l in range(args.start, args.stop + 1):
        cseq = seq[0:l]
        cname = "{}-len_{}".format(_sname, l)

        csfile = ril.sys_suboptimals(cname, cseq,
                                     RNAsubopt=args.RNAsubopt,
                                     ener=args.s_ener,
                                     temp=args.temperature,
                                     noLP=args.noLP,
                                     sort=[
                                         '|', 'sort', '-T', args.s_sortdir, '-k3r', '-k2n'],
                                     verb=args.verbose,
                                     force=args.force)
        sfiles.append(csfile)
        prog.inc()
    return sfiles


def set_p0(bfile, l, lastlines, curlmin, newlmin, cutoff, verb):
    """
    remap densities from :lastlines: using the information in plist
    """

    lminmap = c.defaultdict(int)
    for x, y in zip(curlmin, newlmin):
        lminmap[x] = y

    p0dict = c.defaultdict(float)
    for (i, pop) in enumerate(lastlines[0].split()):
        if i == 0:
            time = float(pop)
        else:
            p0dict[lminmap[i]] += float(pop)
            if float(pop) > cutoff:
                ss, en = get_structure(bfile, i, energy=True)
                print "{:3d} {:3d} {:f} {:s} {:6.2f} {:4d} => {:d}".format(
                    l, i, float(pop), ss, float(en), i, lminmap[i])
                if lminmap[i] == 0:
                    raise LostPopulationError('Lost significant population!')

    p0 = []
    p0sum = 0.0
    for (x, y) in p0dict.items():
        if x == 0:
            continue
        if y > cutoff:
            p0.append("{:d}={:f}".format(x, y))
            p0sum += y

    if verb:
        print "# Total population {:.3f}\n".format(p0sum)
    return p0


def get_structure(bfile, idx, energy=False):
    ss = ''
    with open(bfile, 'r') as bar:
        for n, line in enumerate(bar):
            if n == idx:
                [ss, en] = line.strip().split()[1:3]
                break
    if energy:
        return ss, en if ss else ''
    else:
        return ss


def get_mapping_dict(oldbar, minfo):
    '''
    * parse the old barfile
        o = [[idx, struct],...]
    * parse the current barfile
        c = [[idx, struct],...]
    * dictionary[old_idx]=[new_idx,dist]
    '''
    mapper = c.defaultdict(int)
    with open(oldbar, 'r') as old, \
            open(minfo, 'r') as info:

        mapinfo = []
        for n, line in enumerate(info):
            if re.match('[\.\(\)]+', line.strip()):
                gstr, sptidx, energy, fmin, fminT, gmin, gminT = line.strip().split()
                mapinfo.append([gminT, gstr])
            elif re.match('not in hash', line.strip()):
                print "# structure not in hash"
                mapinfo.append([0, ''])
            elif re.match('not yet assigned', line.strip()):
                print "# structure not yet assigned"
                mapinfo.append([0, ''])

        for n, line in enumerate(old):
            if n == 0:
                continue
            else:  # old idx = enumerate n
                [old_idx, old_min] = line.strip().split()[0:2]
                if n > len(mapinfo):
                    raise Exception('ERROR: To many lines in barfile: '
                                    + minfo + ' missing mapping information!')
                mapper[int(old_idx)] = int(mapinfo[n - 1][0])

    return mapper


def pathlist(mapdata):
    '''
    Translate a list of dictionaries with mapping information into
    a list of trajectories for lmin-transitions

    :return: plist
    '''

    # initialize path-lists with first lmins
    plist = [[i + 1] for i in range(len(mapdata[0]))]

    # append the mapping for each current minimum
    for m, mdict in enumerate(mapdata):
        seen = set()
        for path in plist:
            last = path[-1]
            # print path, "last: ", last, mdict[last]
            path.append(mdict[last])
            seen.add(mdict[last])

        # for path in plist: print "update: ", path

        # add new pathlists for lmins that appear in
        # the next landscape for the first time
        if m == len(mapdata) - 1:
            break
        for newmin in sorted(
                mapdata[m + 1], key=mapdata[m + 1].get, reverse=False):
            if newmin not in seen:
                newpath = [0 for i in range(len(plist[0]) - 1)]
                newpath.append(newmin)
                # print "Adding new path: ", newpath
                plist.append((newpath))

    return sorted(plist, key=lambda m: m[-1])


def add_barmap_args(parser):
    """ A collection of arguments that are used by BarMap """
    ril.argparse_add_arguments(parser,
                               RNAsubopt=True,
                               barriers=True,
                               treekin=True,
                               noLP=True, temperature=True,
                               tmpdir=True, name=True, force=True, verbose=True,
                               start=True, stop=True, k0=True, tX=True, cutoff=True)

    parser.add_argument("--plot_title", default='')

    parser.add_argument("--pyplot", action="store_true",
                        help="Plot the simulation using matplotlib. Interpret the legend \
          using the *log* output")
    parser.add_argument("--xmgrace", action="store_true",
                        help="Print a plot for xmgrace. " +
                        "Interpret the legend using the *log* output")

    parser.add_argument("--adaptive", action="store_true",
                        help="Automatically raise suboptimal energy range if computations fail.")

    parser.add_argument("--s_sortdir", default="/tmp", action="store",
                        help=argparse.SUPPRESS)
    return


def main(args):
    """ BarMap-v2.0 -- cotransriptional folding
      Dependencies: RNAsubopt, barriers, treekin

      The implementation is split into 4 steps ...
    """
    # Read Input & Update Arguments
    name, seq = ril.parse_vienna_stdin(sys.stdin)

    # One name, just to be clear ...
    (args.name, name) = (args.name, args.name) if args.name else (name, name)

    if args.stop is None:
        args.stop = len(seq)
    else:
        seq = seq[:args.stop]

    print "# Input: {:s} {:s}".format(name, seq)

    if args.s_ener is None:
        args.s_ener, args.s_maxn = ril.sys_subopt_range(seq,
                                                        nos=args.s_maxn, maxe=args.s_maxe, verb=(args.verbose > 0))
        print "# Energyrange {:.2f} computes {:d} sequences".format(
            args.s_ener, args.s_maxn)
    elif args.verbose:
        args.s_ener, args.s_maxn = ril.sys_subopt_range(seq,
                                                        nos=0, maxe=args.s_ener, verb=False)
        print "# Energyrange {:.2f} computes {:d} sequences".format(
            args.s_ener, args.s_maxn)

    if not args.tmpdir:
        args.tmpdir = 'BarMap_' + args.name

    if not os.path.exists(args.tmpdir):
        os.makedirs(args.tmpdir)

    """# Starting with BarMap computations ... """

    while True:
        print """# writing RNAsubopt files ... """
        sname = "{}/{}-ener_{:.2f}".format(args.tmpdir, args.name, args.s_ener)
        #if args.circ: myfile += '_circ'
        if args.noLP:
            sname += '_noLP'
        sfiles = barmap_subopts(sname, seq, args)

        print """# writing barriers files ... """
        bname = "{}-minh_{}-maxn_{}-k0_{}".format(sname,
                                                  args.b_minh, args.b_maxn, args.k0)
        bfiles = barmap_barriers(bname, seq, sfiles, args)

        print """# writing/parsing mapping information ... """
        plist = barmap_mapping(bname, seq, args)

        print """# simulations using treekin ... """
        try:
            tfiles = barmap_treekin(bname, seq, bfiles, plist, args)
            break
        except LostPopulationError as e:
            if args.adaptive:
                args.s_ener += 2
                print Warning('repeating caluclations with higher energy:', args.s_ener)
            else:
                print Warning('caluclations failed with current suboptimal energy range:', args.s_ener)
                break
        except SubprocessError as e:
            if args.adaptive:
                args.s_ener += 2
                print Warning('repeating caluclations with higher energy:', args.s_ener)
            else:
                print Warning('caluclations failed with current suboptimal energy range:', args.s_ener)
                break

    if args.xmgrace or args.pyplot:
        print """# Processing treekin results for plotting ... """
        courses = get_plot_data(tfiles, plist, args)

        if args.xmgrace:
            grfile = plot_xmgrace(courses, plist, args)
            print "# Your results have been plotted in the file: {}".format(grfile)

        if args.pyplot:
            plotfile = plot_matplotlib(name, seq, courses, plist, args)
            print "# Your results have been plotted in the file: {}".format(plotfile)

    print "# Thank you for using BarMap b(^.^)d"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        # formatter_class=argparse.RawTextHelpFormatter,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # formatter_class=argparse.MetavarTypeHelpFormatter,
        description='echo sequence | %(prog)s [options]')

    add_barmap_args(parser)

    args = parser.parse_args()

    main(args)
