# ribolands 
energy landscapes and folding kinetics of nucleic acids

**ribolands** is a python package to compute folding kinetics on dynamic or
bimolecular energy landscapes. First of all, it provides an interface to the
programs `RNAsubopt`, `barriers` and `treekin`, which have to installed
separately. See below for an example workflow.

Two scripts for cotranscriptional folding are part of the `ribolands` package: 

  * **BarMap**: folding kinetics on dynamic energy landscapes. For each
    sequence length, the coarse-grained barriers landscape is computed. During
    kinetic simulations, a mapping between subsequent landscapes is used to
    transfer occupancy from one landscape to the next. This is mostly a
    reimplementation of `BarMap` by [Hofacker et al. (2010)].
    ```sh
    echo "CUGCGGCUUUGGCUCUAGCC" | BarMap.py -v --plot_linlog 
    ```

  * **DrTransformer**: Short for "DNA-to-RNA Transformer", the program
    computes cotranscriptional folding of larger RNAs by generating a
    conformation graph of local-minimum conformations at every transcription
    step. Importantly, only conformations reachable from the previous
    conformation graph are included, rates between these conformations are
    computed using the `findpath` direct path heuristic, see [Flamm et al.
    (2001)]. `DrTransformer` uses `treekin` to calculate folding kinetics.
    ```sh
    echo "CUGCGGCUUUGGCUCUAGCC" | DrTransformer.py -v --pyplot
    ```

## ViennaRNA dependencies
`ribolands` uses [RNAsubopt] from the [ViennaRNA package], [barriers] and
[treekin] for landscape computations. Make sure that you have the latest
versions installed, i.e. `treekin-v0.4`, `barriers-v1.6` and, recommended,
`ViennaRNA-v2.2` or later.

### Examples
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
  - Stefan Badelt. [Control of RNA function by conformational design.] PhD thesis, University of Vienna, 2016
  - BarMap: RNA folding on dynamic energy landscapes [Hofacker et al. (2010)] 

## Installation
```sh
  python setup.py install
```

### local installation
```sh
  python setup.py install --prefix=/your/destination/ribolands-x.y.r
```
Do not forget to set your environment variables when using local installations:
```sh
  export PATH=/your/destination/ribolands-x.y.r/bin:$PATH
  export PYTHONPATH=/your/destination/ribolands-x.y.r/lib/python2.7/site-packages/:$PYTHONPATH
```
  
## Version
0.2.0

## Development / Unittests
  python setup.py test

## Build the documentation
  sphinx-build -b html docs ~/your/html/sourcedir/

## License
MIT

[//]: References
[Hofacker et al. (2010)]: <http://dx.doi.org/10.1261%2Frna.2093310>
[Flamm et al. (2001)]: <http://rnajournal.cshlp.org/content/7/2/254.short>

[ViennaRNA package]: <http://www.tbi.univie.ac.at/RNA>
[RNAsubopt]: <http://www.tbi.univie.ac.at/RNA/RNAsubopt.1.html>
[barriers]: <http://www.tbi.univie.ac.at/RNA/Barriers>
[treekin]: <http://www.tbi.univie.ac.at/RNA/Treekin>

[Control of RNA function by conformational design.]: <http://www.tbi.univie.ac.at/newpapers/pdfs/TBI-t-2016-1.pdf>

[ribolands]: <https://www.tbi.univie.ac.at/RNA/ribolands>

