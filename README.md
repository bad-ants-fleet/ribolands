#### 'ribolands': folding kinetics using ViennaRNA ####

Full documentation should at some point be available at: [http://www.tbi.univie.ac.at/RNA/ribolands/](http://www.tbi.univie.ac.at/RNA/ribolands/)

#### Build the documentation ####

sphinx-build -b html docs ~/your/html/sourcedir/

##### Installation #####

python setup.py install

#### local installation ####
python setup.py install --prefix=/your/destination/ribolands-x.y.r

**Do not forget to set your environment variables when using local installations**

export PATH=/your/destination/ribolands-x.y.r/bin:$PATH
export PYTHONPATH=/your/destination/ribolands-x.y.r/lib/python2.7/site-packages/:$PYTHONPATH


