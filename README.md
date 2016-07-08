# ribolands -- energy landscapes and folding kinetics of nucleic acids

**ribolands** is a python interface to ViennaRNA style folding kinetics. [...]
Full documentation will be available at []

### Examples
```sh
>>> from ribolands import sys_barriers

>>> files = sys_barriers('ACGUGACUG', maxn=50, minh=1.0)

```

## Installation
  python setup.py install

### local installation
  python setup.py install --prefix=/your/destination/nuskell-x.y.r
  
Do not forget to set your environment variables when using local installations:
  
  export PATH=/your/destination/nuskell-x.y.r/bin:$PATH
  export PYTHONPATH=/your/destination/nuskell-x.y.r/lib/python2.7/site-packages/:$PYTHONPATH
  
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

