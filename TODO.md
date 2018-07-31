# TODO list

### General style updates:
  - fix preabmles of all files (unicode? utf-8?)
  - check dependencies
  - print verbose output to STDOUT (instead of STDERR)
  - write unittests
  - compile manpage
  - write changelog

### ribolands.syswraps
  - encapsulate syswraps into tmp-directory (bec of barriers)

### other
  - replace spatch.py with ViennaRNA 2.2 soft-constraints
  - add scripts/cofold_barriers.py
  - write a decorator function to submit sys_jobs to the SGE, 
      rather than computing them locally

### BarMap:
Calculate BarMap sequentially, take all minima that were occupied in the last
round, check the saddle points and then use the energy of the saddle point to
determine the subopt range, warn if subopt range exceeds some threshold.



