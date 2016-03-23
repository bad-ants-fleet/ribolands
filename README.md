#### 'vrna.pylands': folding kinetics using ViennaRNA ####

Full documentation should at some point be available at: [http://www.tbi.univie.ac.at/RNA/pylands/](http://www.tbi.univie.ac.at/RNA/pylands/)

#### Build the documentation ####

sphinx-build -b html doc ~/your/html/sourcedir/

##### Installation #####

python setup.py install

#### local installation ####
python setup.py install --prefix=/your/destination/rnaworld-x.y.r

**Do not forget to set your environment variables when using local installations**

export PATH=/your/destination/rnaworld-x.y.r/bin:$PATH
export PYTHONPATH=/your/destination/rnaworld-x.y.r/lib/python2.7/site-packages/:$PYTHONPATH


