# ribolands
energy landscapes and folding kinetics of nucleic acids

**ribolands** is a python interface to compute folding kinetics on dynamic or
bimolecular energy landscapes. First of all, it provides a python interface to
the programs `RNAsubopt`, `barriers` and `treekin`, which have to installed
separately. See below for an example workflow.

Two scripts for cotranscriptional folding are part of the `ribolands` package: 

  * **BarMap.py**: folding kinetics on dynamic energy landscapes. For each
    sequence length, the coarse-grained barriers landscape is computed. During
    kinetic simulations, a mapping between subsequent landscapes is used to
    transfer occupancy from one landscape to the next.

  * **DrTransformer.py**: cotranscriptional folding on a reduced conformation
    graph.

### Examples
```sh
>>> from ribolands import sys_subopt_range, sys_suboptimals, sys_barriers

>>> [name, seq] = ['test', 'ACUGAGGUCGAU']

>>> ener, nos = sys_subopt_range(seq, nos=100000)
>>> sfile = sys_suboptimals(name, seq, ener=ener)
>>> [sfile, bfile, efile, rfile, psfile] = sys_barriers(name, seq, sfile, maxn=50, minh=1.0, rates=True)
>>> tfile = sys_treekin(name, seq, bfile, rfile, p0=['2=1'], t0=1e-6, t8=1e10)
```

## Installation
  python setup.py install

### local installation
  python setup.py install --prefix=/your/destination/ribolands-x.y.r
  
Do not forget to set your environment variables when using local installations:
  
  export PATH=/your/destination/ribolands-x.y.r/bin:$PATH

  export PYTHONPATH=/your/destination/ribolands-x.y.r/lib/python2.7/site-packages/:$PYTHONPATH
  
## Version
0.2.0

### Devolment / Unittests
  python setup.py test

### Build the documentation
  sphinx-build -b html docs ~/your/html/sourcedir/

### Todos

 - Write Tests
 - Sphinx Documentation
 - Release

### License
MIT

[//]: References
[ribolands]: <https://rna.tbi.univie.ac.at/ribolands>
