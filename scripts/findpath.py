#!/usr/bin/env python3

# Written by Stefan Badelt (stef@tbi.univie.ac.at)

import sys
import argparse
import RNA

def costruct(seq, cut = None):
    """ Translate between 'input' and 'internal' RNAcofold sequences """
    delim = '&'
    if delim in seq :
        cut = seq.index(delim) + 1
        table = str.maketrans(dict.fromkeys('&'))
        seq = seq.translate(table)
    elif cut and cut != -1 :
        seq = seq[:cut-1]+delim+seq[cut-1:]
        cut = None
    return seq, cut

def findpath_wrap(fc, s1, s2, maxE, fpath, cut = None, verbose=False):
    """A wrapper for ViennaRNA findpath functions.

    Args:
        fc (RNA.fold_compound()): 
    
    """
    sE = None
    if verbose :
        if maxE:
            dcal_bound = int(round(maxE * 100))
            path = fc.path_findpath(s1, s2, maxE=dcal_bound, width=fpath)
        else:
            path = fc.path_findpath(s1, s2, width=fpath)

        if len(path):
            sE = round(fc.eval_structure(s1), 2)
            for e, step in enumerate(path):
                ss, _ = costruct(step.s, cut)
                print("{:2d} {:s} {:6.2f}".format(e, ss, step.en))
                if step.en > sE: sE = step.en
    else :
        if maxE:
            dcal_bound = int(round(maxE * 100))
            dcal_sE = fc.path_findpath_saddle(s1, s2, maxE=dcal_bound, width=fpath)
        else :
            dcal_sE = fc.path_findpath_saddle(s1, s2, width=fpath)
        if dcal_sE is not None:
            sE = float(dcal_sE)/100
    return sE

def main():
    """A wrapper for the ViennaRNA findpath routine."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-v","--verbose", action="store_true",
            help="Verbose output: print the folding pathway.")
    parser.add_argument("-w","--search-width", type=int, default=100,
            help="Adjust upper bound for findpath search.")
    parser.add_argument("-m","--max-energy", type=float, default=None,
            help="Specify upper bound for barrier energy (kcal/mol).")
    args = parser.parse_args()

    for e, line in enumerate(sys.stdin):
        if e == 0 :
            seq = line.strip()
        elif e == 1 :
            ss1, cut = costruct(line.strip())
        elif e == 2 :
            ss2, cut_ = costruct(line.strip())

    assert cut == cut_

    md = RNA.md()
    fc = RNA.fold_compound(seq, md)

    if args.verbose:
        print("   {}".format(seq))

    saddle = findpath_wrap(fc, ss1, ss2, args.max_energy, 
            args.search_width, cut, args.verbose)

    if saddle is not None:
        e1 = round(fc.eval_structure(ss1), 2)
        barrier = saddle-e1
        print("Saddle: {:6.2f} kcal/mol | Barrier: {:6.2f} kcal/mol | Search width paramter: {}".format(saddle, barrier, args.search_width))
    else:
        print("No valid saddle found below energy: {:6.2f} kcal/mol | Search width parameter: {}".format(args.max_energy, args.search_width))


if __name__ == '__main__':
  main()

