# ribolands - kinetic analysis of nucleic acids

**ribolands** is a Python wrapper for multiple nucleic acid kinetic analysis
programs, currently `Kinfold`, `RNAsubopt`, `barriers` and `treekin`. Those 
programs have to installed separately. See below for an example workflow.

Two scripts for cotranscriptional folding analysis are installed by `ribolands`: 

  * **BarMap**: Folding kinetics on dynamic `barriers`-type energy landscapes.
    For each sequence length, suboptimal structures are enumerated, and the
    landscape is coarse-grained into basins of a minimal height (using the
    program barriers). For kinetic simulations, a mapping between subsequent
    landscapes is generated to transfer occupancy from one landscape to the
    next.  This is mostly a reimplementation of `BarMap` by [Hofacker et al.
    (2010)], but it makes use of more recent functionality of `barriers` and
    `treekin`.
    ```sh
    echo "CUGCGGCUUUGGCUCUAGCC" | BarMap.py --pyplot
    ```

  * **DrTransformer**: Short for "DNA-to-RNA Transformer". A program for
    cotranscriptional folding simulations of larger RNAs, using a heuristic to
    approximate the relevant structures at every transcription step.  
    `DrTransformer` uses the `ViennaRNA package` to calculate transition rates
    between relevant structures and `treekin` to simulate folding kinetics.
    ```sh
    echo "CUGCGGCUUUGGCUCUAGCC" | DrTransformer.py
    ```

## Installation
```sh
  ~$ python setup.py install
  ~$ python -m pytest tests/ -v -s
```

## ViennaRNA dependencies
ribolands uses [RNAsubopt] from the [ViennaRNA package], as well as [barriers]
and [treekin] for landscape computations. Make sure that you have recent
versions installed: 
 - `ViennaRNA>=v2.4.13`
 - `treekin>=v0.5.1`
 - `barriers>=v1.8.1` 

### Examples (TODO, outdated)
```
>>> from ribolands import sys_subopt_range, sys_suboptimals, sys_barriers, sys_treekin

>>> [name, seq] = ['test', 'ACUGAGGUCGAU']

>>> ener, nos = sys_subopt_range(seq, nos=100000)
>>> sfile = sys_suboptimals(name, seq, ener=ener)
>>> [sfile, bfile, efile, rfile, psfile] = sys_barriers(name, seq, sfile, maxn=50, minh=1.0, rates=True)
>>> [tfile, efile] = sys_treekin(name, seq, bfile, rfile, p0=['2=1'], t0=1e-6, t8=1e10)
```

## Cite
If you are using `BarMap` or `DrTransformer` please cite: 
  - Stefan Badelt, ... in prepration. 
  - BarMap: RNA folding on dynamic energy landscapes [Hofacker et al. (2010)] 
 
## Version
0.9

## License
MIT

[//]: References
[Hofacker et al. (2010)]: <http://dx.doi.org/10.1261%2Frna.2093310>
[Flamm et al. (2001)]: <http://rnajournal.cshlp.org/content/7/2/254.short>
[ViennaRNA package]: <http://www.tbi.univie.ac.at/RNA>
[RNAsubopt]: <http://www.tbi.univie.ac.at/RNA/RNAsubopt.1.html>
[barriers]: <http://www.tbi.univie.ac.at/RNA/Barriers>
[treekin]: <http://www.tbi.univie.ac.at/RNA/Treekin>
[ribolands]: <https://www.tbi.univie.ac.at/RNA/ribolands>

