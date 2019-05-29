# Ribolands - Todo list

### General style updates
  - fix preabmles of all files (utf-8?)
  - print verbose output to STDOUT (instead of STDERR)
  - write unittests
  - compile manpage
  - write changelog

### Checklist
  - check dependencies
  - ./setup.py test
  - pytest tests/

### ribolands.syswraps
  - encapsulate syswraps into tmp-directory (bec of barriers)
  - consider updating to pourRNA instead of barriers
  - use yield / next to enforce structure on subopt / barriers / treekin

### scripts
  - DrTransformer.py v1.0 (finalize)
    - DrProfile.py v1.0     (finalize)
  - BarMap.py v1.0        (finalize)
  - spatch.py v0.1        (remove or update to ViennaRNA soft constraints)
  - cobarriers.py v0.1    (finalize or remove)

### BarMap improvements:
Calculate BarMap sequentially, take all minima that were occupied in the last
round, check the saddle points and then use the energy of the saddle point to
determine the subopt range, warn if subopt range exceeds some threshold.


# Repository structure:
Should DrTransformer be part of ribolands? Or should it be dependent on ribolands?

yes, it should be part of ribolands:
    ribolands without DrTransformer means also no BarMap ... but that is sad.
    Why would one maintain ribolands if the most important applications are gone?
    If ribolands is not maintained, DrTransformer will soon stop working...
    => it cannot *just* depend on ribolands, it has to duplicate ribolands

    DrTransformer remains a Python script, that connects treekin, ViennaRNA,
    and networkx. It could be viewed as a pretty good prototype (just like
    BarMap), and so long it is fine to be part of the ribolands project.

    the extended feature library of ribolands make it more flexible to build new
    hybrid approaches. E.g. KinKong, or BarMap/Kinfold/ stuff... while giving
    different perspectives on problems, improvements of the code base will help
    all of these programs.

no, it should be its own program:
    DrTransformer will be published and it is the main project, so why would
    you hide it in a package with a different name?

    DrTransformer is more difficult to maintain if you always have to be
    careful that you do not destroy other projects like BarMap... 
    (although... that's what unittests are for, right?)

Should Ribolands be hosted at ViennaRNA GitHub?
    Yes, there should be two types of dependencies: ViennaRNA-group programs,
    and pip installable python projects. 
    

