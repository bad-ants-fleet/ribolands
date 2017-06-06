Cotranscriptional folding using the ribolands package
=====================================================

DrTransformer
-------------
Short for "DNA-to-RNA Transformer", the program computes cotranscriptional
folding of larger RNAs by generating a conformation graph of local-minimum
conformations at every transcription step. Importantly, only conformations
reachable from the previous conformation graph are included, rates between
these conformations are computed using the `findpath` direct path heuristic,
see [Flamm et al. (2001)]. `DrTransformer` uses `treekin` to calculate folding
kinetics.

  ``$ echo "CUGCGGCUUUGGCUCUAGCC" | DrTransformer.py --pyplot``



BarMap
------
`BarMap` computes folding kinetics on dynamic energy landscapes. For each
sequence length, the coarse-grained barriers landscape is computed. During
kinetic simulations, a mapping between subsequent landscapes is used to
transfer occupancy from one landscape to the next. This is mostly a
reimplementation of `BarMap` by [Hofacker et al. (2010)].

  ``$ echo "CUGCGGCUUUGGCUCUAGCC" | BarMap.py --pyplot``



