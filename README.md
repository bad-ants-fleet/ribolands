# ribolands - kinetic analysis of nucleic acids

**ribolands** is a Python wrapper for multiple nucleic acid kinetic analysis
programs, currently `Kinfold`, `RNAsubopt`, `barriers` and `treekin`. Those 
programs have to installed separately. 

:warning: Dear users of the cotranscriptional folding simulator
**DrTransformer**, this package provided only a **preliminary version** of the
software which now has been removed. The newest version of **DrTransformer** is
not publicly available yet, but feel free to reach out for beta-testing.
Also, the executable for the cotranscriptional folding software `BarMap` has
been removed (since it is has never been officially released). `BarMap` can be
coded with minimal effort using the ribolands package, feel free to reach out
if you need help.

## Dependencies
`ribolands` uses [RNAsubopt] from the [ViennaRNA package], as well as [barriers]
and [treekin] for landscape computations. Make sure that you have recent
versions installed: 
 - `ViennaRNA>=v2.4.13`
 - `treekin>=v0.5.1`
 - `barriers>=v1.8.1` 


## Installation
```sh
  ~$ python setup.py install
  ~$ python -m pytest tests/ -v -s
```

## Version
0.10 -- removed executables
 - removed DrTransformer, which will be available soon as a separate python package.
 - removed all other executables in order to focus on library development.
 - lot's of interface changes to build workflows (work in progress)

## License
MIT

[//]: References
[Hofacker et al. (2010)]: <http://dx.doi.org/10.1261%2Frna.2093310>
[ViennaRNA package]: <http://www.tbi.univie.ac.at/RNA>
[RNAsubopt]: <http://www.tbi.univie.ac.at/RNA/RNAsubopt.1.html>
[barriers]: <http://www.tbi.univie.ac.at/RNA/Barriers>
[treekin]: <http://www.tbi.univie.ac.at/RNA/Treekin>
[ribolands]: <https://www.tbi.univie.ac.at/RNA/ribolands>
[drtransformer]: <https://www.github.com/ViennaRNA/drtranformer> 

