RNA folding kinetics tutorial using ribolands
=============================================

This section provides a short and concise example-workflow of how to use `ribolands`.
It describes how to (i) compute coarse-grained RNA energy landscapes and (ii)
compute folding kinetics on these landscapes.  It is highly recommended that
users find their own approach by looking deeper into the options of the
different functions. For algorithms and assumptions behind the ``RNAsubopt |
barriers | treekin`` pipeline, please study the respective manpages and the
publications.

workflow for single RNA folding kinetics
----------------------------------------

Import and initialization

>>> from ribolands import sys_subopt_range, sys_suboptimals, sys_barriers, sys_treekin 
>>> [name, seq]  = ["xbix", "CUGCGGCUUUGGCUCUAGCC"]

First, we need to determine an appropriate energy range for suboptimal
structure computation. In some cases, e.g. for short molecules, you might have
a particular energy range that you are satisfied with, but more often it is a
trade-off between the number of secondary structures and the computation time
of ``RNAsubopt`` and ``barriers``. Hence, in practice it is useful to
precompute an energy range that results in a certain number of seconary
structures. As a rule of thumb, 500.000 will lead to realtively quick
results, while 10.000.000 can take a few hours of computation time. For
larger values, it is recommended to compile ``barriers`` using non-standard
options.
  
>>> ener, nos = sys_subopt_range(seq, nos=500000)
>>> print "Energy range: {%.2f} computes {%d} sequences".format(ener, nos)

The precomputed energy range lets us savely call ``RNAsubopt``, without bad
surprises during subsequent calculations in the pipeline. Many of the standard
``RNAsubopt`` commandline options are available as optional arguments [source].
Note, that not all of the parameters have the same defaults values, for example,
the output is `sorted` by default, in order to be compatible with ``barriers``.
Also, uncompressed subopt-results can easily grow to a few gigabites;
by default ``ribolands`` operates on subopt files compressed using ``gzip``.

>>> sfile = sys_suboptimals(name, seq, ener=ener)
>>> print "Results from RNAsubopt have been printed in the file:", sfile

.. try::
.. 
..   $ zcat xbix.spt.gz | head

Next, the suboptimal structures are coarse-grained into an energy landscape.
Again, standard parameters in ``ribolands`` differ from those in the `barriers`
excutable, e.g. we use `single-base-pair moves` as the default. In the
following, we compute landscapes with a minimum basin-height of 1.0 kcal/mol and
terminate after a maximum of 50 local minimum conformations. The option to
compute rates between macro-states is switched on, as we need them for
subsequent ``treekin`` computations.

>>> [sfile, bfile, efile, rfile, psfile] = sys_barriers(name, seq, sfile, maxn=50, minh=1.0, rates=True)
>>> print "Standard output of barriers printed to:", bfile
>>> print "Standard error of barriers printed to:", efile
>>> print "Rate matrix printed to:", rfile
>>> print "Projection of the energy landscape into a barrier-tree:", psfile

To compute folding kinetics on the coarse-grained energy landscape, you need to
pass both the `bfile` and the `rfile` on to `treekin`. Use the option p0 to
specify one or more starting conformations and their occupancy. With t0 and t8
we set the start and end time of the simulation in `arbitrary time units`.

>>> tfile = sys_treekin(name, seq, bfile, rfile, p0=['2=1'], t0=1e-6, t8=1e10)
>>> print "A file with RNA folding trajectories", tfile

Use a standard plotting tool to display the results in the treekin output file.


