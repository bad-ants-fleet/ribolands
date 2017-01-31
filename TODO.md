# Issues

### General checklist:
  - fix preabmles of all files
  - check dependencies
  - print verbose output to STDOUT or LOGFILE (instead of STDERR)
  - write unittests
  - compile manpage

### ribolands.syswraps
  - encapsulate syswraps into tmp-directory (bec of barriers)
  - find a nicer way to deal with k0 (e.g. treekin --times)

### other
  - replace spatch.py with ViennaRNA 2.2 soft-constraints
  - add crnsimulator
  - add scripts/cofold_barriers.py
  - write a decorator function to submit sys_jobs to the SGE, 
      rather than computing them locally

### Projects
  - finish the binding pocket implementations (XOR, Lucs)
  - show a folding path with a partiuclar Luck's visualization?


