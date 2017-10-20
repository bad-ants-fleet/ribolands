
import re
import math
import networkx as nx
import subprocess as s
from struct import pack

import RNA
import ribolands as ril

class TrafoUsageError(Exception):
  pass

class TrafoAlgoError(Exception):
  pass

class DebuggingAlert(Exception):
  pass

class ConformationGraph(nx.DiGraph):

  def __init__(self, fullseq, vrna_md, RT = 0.61632077549999997):
    """ Conformation Graph initialization.
    """
    super(ConformationGraph, self).__init__()

    self._full_sequence = fullseq
    self._model_details = vrna_md
    self._fold_compound = RNA.fold_compound(fullseq, vrna_md)
    self._RT = RT

    # Private instance variables:
    self._transcript_length = 0
    self._total_time = 0
    self._nodeid = 0

    # Default parameters:
    self._cutoff = 0.01
    self._fpath = 20
    self._maxdG = None
    self._k0 = 2e5

    # expand() 
    #self._mfree = 6 #args.min_breathing
    #self._exp_mode = 'default' 

    # Output formats:
    # TODO: plot only a certain region "start - end parameter"
    self._matplot = None 
    self._xmgrace = None 
    self._drforna = None 
    self._logfile = None 

  @property
  def full_sequence(self):
    return self._full_sequence

  @property
  def transcript(self):
    return self._full_sequence[0:self._transcript_length]

  @property
  def occupancy_cutoff(self):
    return self._cutoff

  @occupancy_cutoff.setter
  def occupancy_cutoff(self, occu):
    self._occupancy_cutoff = occu

  @property
  def findpath_width(self):
    return self._fpath

  @findpath_width.setter
  def findpath_width(self, val):
    self._fpath = val

  @property
  def logfile(self):
    return self._logfile

  @logfile.setter
  def logfile(self, fh):
    """ filehandle or None to set a logfile """
    self._logfile = fh

  def get_saddle(self, s1, s2):
    if self.has_edge(s1,s2) :
      return self[s1][s2]['saddle']
    else :
      return None

  def get_rate(self, s1, s2):
    if self.has_edge(s1,s2) :
      return self[s1][s2]['weight']
    else :
      return None

  def add_transition_edge(self, s1, s2, ts = None):
    """ compute the transition rates between two structures: s1 <-> s2, 
      where s2 is always the new, energetically better structure.
      NOTE: actually, that is not always the case if you use exterior_folding,
      then, s1 is always the conformation that is already present in the confomration
      graph...
    """
    maxdG = self._maxdG
    fpath = self._fpath
    _RT = self._RT
    k0 = self._k0

    fullseq = self._full_sequence
    md = self._model_details
    fc = self._fold_compound
  
    # Minimum between direct and in-direct path barriers
    if ts : # then we know that the indirect path has to be part of saddles
      tsE1 = self.get_saddle(s1,ts)
      tsE2 = self.get_saddle(s2,ts)
      tsE = max(tsE1, tsE2)

    saddleE = self.get_saddle(s1,s2)
    if saddleE is None :
      if ts:
        sE = float(fc.find_path_once(s1, s2, int(round(tsE * 100)), fpath))/100
        saddleE = min(tsE, sE)
        sE = float(fc.find_path_once(s2, s1, int(round(saddleE * 100)), fpath))/100
        saddleE = min(saddleE, sE)
      else :
        saddleE = float(fc.path_findpath_saddle(s1, s2, fpath))/100
  
    # invalid neighbor saddle energies are not stored for future calls.
    if maxdG :
      valid = (ts is not None or saddleE-self.node[s1]['energy'] <= maxdG)
    else :
      valid = True
  
    if valid :
      e1 = self.node[s1]['energy']
      e2 = round(fc.eval_structure(s2), 2)
  
      saddleE = max(saddleE, max(e1,e2)) # ensure saddle is not lower than s1, s2
  
      # Energy barrier
      dG_1s = saddleE-e1
      dG_2s = saddleE-e2
  
      # Metropolis Rule
      k_12 = k0 * math.exp(-dG_1s/_RT)
      k_21 = k0 * math.exp(-dG_2s/_RT)
  
      self.add_weighted_edges_from([(s1, s2, k_12)])
      self.add_weighted_edges_from([(s2, s1, k_21)])
      self[s1][s2]['saddle'] = saddleE
      self[s2][s1]['saddle'] = saddleE
  
      #print "#Added Edge:", s1, s2, "({}, {:g}, {:g})".format(valid, k_12, k_21)
 
    return valid

  def expand(self, extend = 1, exp_mode = 'default', mfree = 6, cutoff = None):
    """ Find new neighbors and add them to the Conformation Graph

    The function is devided into two parts. 1) The current mfe structure
    is connected to all present structures, 2) The conformation graph is
    expanded using helix-breathing.

    Args:
      extend (int, optional): Number of nucleotide extensions. Defaults to 1.
      exp_mode (str, optional): Expansion mode: choose from "mfe-only": only
        use current mfe as potential new neighbor. "breathing-only": only use
        breathing neighborhood. "default": do both mfe and breathing.
      mfree (int, optional): minimum number of freed bases during a
        helix-opening step. Defaults to 6.
      cutoff (flt, optional): Occupancy cutoff for neighbor generation.
        Defaults to None, and therefore the global cutoff set in the
        ConformationGraph 

    # TODO: coarse grain results... that means if a rate is much faster than
    # 1/t8 it will equilibrate, so lets just distribute the occupancy there upfront.
    # This reduces the output for visualization as well...

    Returns:
      int: Number of new nodes
    """
    fseq = self._full_sequence

    if cutoff is None:
      cutoff = self._cutoff

    self._transcript_length += extend
    if self._transcript_length > len(fseq):
      self._transcript_length = len(fseq)
    seq = self.transcript

    csid = self._nodeid
    md = self._model_details
    fc_full = self._fold_compound

    if exp_mode not in ['default', 'mfe-only', 'breathing-only']:
      raise TrafoUsageError('unknown expansion mode')

    # Add MFE
    fc_tmp = RNA.fold_compound(seq, md)
    ss, mfe = fc_tmp.mfe()
    future = '.' * (len(fseq)-len(seq))
    ss = ss + future
    #print "{}\n{} {:6.2f}".format(seq, ss, mfe)

    # If there is no node bec we are in the beginning, add the node,
    # otherwise, go through all nodes and try to add transition edges
    if len(self) == 0 :
      en = round(fc_full.eval_structure(ss), 2)
      self.add_node(ss, energy=en, occupancy=1.0, 
          identity=self._nodeid, active=True, last_seen=0)
      self._nodeid += 1

    elif exp_mode == 'default' or exp_mode == 'mfe-only':
      # Try to connect MFE to every existing state
      for ni in self.nodes() :
        if self.node[ni]['active'] == False : continue
        if ni == ss : continue
        if self.has_edge(ni,ss) :
          self.node[ss]['active'] = True # in case it was there but inactive
          continue

        if self.has_node(ss): # from a previous iteration
          if self.add_transition_edge(ni, ss):
            self.node[ss]['active'] = True # in case it was there but inactive

        elif self.add_transition_edge(ni, ss):
          en = round(fc_full.eval_structure(ss), 2)
          self.node[ss]['active'] = True
          self.node[ss]['last_seen'] = 0
          self.node[ss]['energy'] = en
          self.node[ss]['occupancy'] = 0.0
          self.node[ss]['identity'] = self._nodeid
          self._nodeid += 1

    if exp_mode == 'default' or exp_mode == 'breathing-only':
      # Do the helix breathing graph expansion

      # Initialize a dictionary to store the feature expansion during
      # each graph expansion round: ext_moves[ext_seq] = [set((con,paren),...), structure]
      # where ext_seq = exterior-loop sequence with ((xxx)) replacing constrained elements
      ext_moves = dict()
      for ni, data in self.nodes_iter(data=True):
        if data['active'] == False : continue
        en  = data['energy']
        occ = data['occupancy']
        if occ < cutoff : continue

        # short secondary structure (without its future)
        sss = ni[0:len(seq)]

        # compute a set of all helix breathing open steps
        opened = open_breathing_helices(seq, sss, mfree)
        #print 'opened', opened

        # do a constrained exterior loop folding for all of them ...
        for onbr in opened :
          nbr, ext_seq = fold_exterior_loop(md, seq, onbr, ext_moves)

          future = '.' * (len(ni) - len(nbr))
          nbr += future

          if ni == nbr: 
            continue

          if self.has_edge(ni, nbr):
            self.node[nbr]['active'] = True # in case it was there but inactive
            continue

          if self.has_node(nbr):
            if self.add_transition_edge(ni, nbr):
              self.node[nbr]['active'] = True # in case it was there but inactive
          elif self.add_transition_edge(ni, nbr) :
            enbr = round(fc_full.eval_structure(nbr), 2)
            self.node[nbr]['energy'] = enbr
            self.node[nbr]['active'] = True
            self.node[nbr]['last_seen'] = 0
            self.node[nbr]['occupancy'] = 0.0
            self.node[nbr]['identity'] = self._nodeid
            self._nodeid += 1
          else :
            print """# WARNING: Could not add transition edge!"""
            # continue? do ext_moves stuff, but don't add to ext_moves?
            # or just force adding?
            raise DebuggingAlert('does this really happen?')

          # And now connect the neighbor with its *historic* transitions:
          # We store the exterior-open neighbor here, that means there are
          # three possible reasons for duplication:
          #   1) the same helix was opened / difference is in historic features
          #   2) different (longer) helix was opened / same historic features
          #   2) different helix / different history
          if ext_moves[ext_seq][0] :
            for (parent, child) in ext_moves[ext_seq][0] :
              assert parent != ni # Parents may never be the same
              # Children can be the same, 
              # if the change is within the helix-breathing-move-set
              # p1 .((((((((.(((.......))).))).)))))......................
              # p2 ..(((((((.(((.......))).))).)))).......................
              # c1 .((((((((.(((.......))).))).)))))......(((......)))....
              # c2 .((((((((.(((.......))).))).)))))......(((......)))....

              if child == nbr : continue
                #raise DebuggingAlert('does this really happen?')

              if self.has_edge(parent, ni) :
                if self.has_node(child) and self.has_node(nbr):
                  if self.has_edge(nbr, child): 
                    continue
                  # TODO: Calculate saddleE from saddleE of parents?
                  # sP = self.get_saddle(ni, parent)
                  # sC1 = round(self.node[child]['energy'] + sP - self.node[parent]['energy'], 2)
                  # sC2 = round(self.node[nbr]['energy'] + sP - self.node[ni]['energy'], 2)

                  if self.add_transition_edge(nbr, child):
                    self.node[nbr]['active'] = True # in case it was there but inactive
                    self.node[child]['active'] = True # in case it was there but inactive
                else :
                  print "should we still find an edge?"
                    
          # Track the final structure, every new identical ext-change will be
          # connected, if the parents were connected.
          ext_moves[ext_seq][0].add((ni, nbr))

    if not self.has_node(ss) or (not self.node[ss]['active']):
      print "# WARNING: ", ss, "[mfe secondary structure not connected]"

    for ni in self.nodes() :
      if self.node[ni]['active'] == False :
        self.node[ni]['last_seen'] += 1
        #print 'last seen:', ni, self.node[ni]['last_seen']
      else :
        if self.node[ni]['last_seen'] > 0:
          print 'revisiting:', ni, self.node[ni]['last_seen']
        self.node[ni]['last_seen'] = 0

      if self.node[ni]['last_seen'] >= 10:
        self.remove_node(ni)
    return self._nodeid-csid

  def simulate(self, t0, t8, tmpfile=None): 
    pass

  def get_simulation_files_tkn(self, name):
    """ Make a barriers + rates output from the current conformation graph.

    The printed files are input files for treekin. *.bar files contain the energy
    barriers to transition between local minima. Although every path eventually
    leads to the MFE structure, it can proceed via an energetically worse
    structure first. This is in contrast to files produced by `barriers`, where
    local minima are always *directly* connected to energetically better local
    minima.   
    """
    seq = self.transcript
    logf = self._logfile

    sorted_nodes = filter(lambda (n,d): d['active'], 
        sorted(self.nodes(data=True), 
          key=lambda x: x[1]['energy'], reverse=False))

    barfile_nodes = filter(lambda d: self.node[d]['active'], 
        sorted(self.nodes(data=False), key=lambda x: self.node[x]['energy'], reverse=False))

    if name :
      bfile = name+'.bar'
      rfile = name+'.rts'
      brfile = rfile+'.bin'
      p0 = []
      
      with open(bfile, 'w') as bar, open(rfile, 'w') as rts, open(brfile, 'w') as brts :
        bar.write("     {}\n".format(seq))
        brts.write(pack("i", len(sorted_nodes)))
        for e, (ni, data) in enumerate(sorted_nodes, 1) :

          # Calculate barrier heights to all other basins.
          nextmin = 0
          barrier = 0
          saddleE = None
          nMsE = set()
          for ee, be in enumerate(barfile_nodes, 1):
            if e == ee :
              continue
            sE = self.get_saddle(be,ni)
            if sE is not None:
              nMsE.add((ee, sE))

          mystr = ' '.join(map(lambda(x,y):'({:3d} {:6.2f})'.format(x,y-data['energy']), 
              sorted(list(nMsE), key=lambda x:x[0])))

          bar.write("{:4d} {} {:6.2f} {}\n".format(e, ni[:len(seq)], data['energy'], 
            mystr))

          if logf is not None:
            line = "{:4d} {:4d} {} {:6.2f} {:6.4f} (ID = {:d})\n".format(
                self._transcript_length, e, ni[:len(seq)], 
                data['energy'], data['occupancy'], data['identity'])
            logf.write(line)

          if data['occupancy'] > 0 :
            p0.append("{}={}".format(e,data['occupancy']))

          trates = []
          rates = []
          for (nj, jdata) in sorted_nodes :
            if self.has_edge(ni,nj) :
              rates.append(self[ni][nj]['weight'])
              trates.append(self[nj][ni]['weight'])
            else :
              rates.append(0)
              trates.append(0)
          line = "".join(map("{:10.4g}".format, rates))
          rts.write("{}\n".format(line))
          for r in trates:
            brts.write(pack("d", r))

    else:
      if logf is not None:
        line = "Distribution of structures at the end:\n"
        line += "          {}\n".format(seq)
        for e, (ni, data) in enumerate(sorted_nodes, 1) :
          line += "LAST {:4d} {} {:6.2f} {:6.4f} (ID = {:d})\n".format(e, 
              ni[:len(seq)], data['energy'], data['occupancy'], data['identity'])
        logf.write(line)
      return 

    return [bfile, brfile, p0, sorted_nodes]

  #   update_time_and_occupancies_tkn(self, tfile)
  def update_occupancies_tkn(self, tfile, sorted_nodes):
    """
      Update the occupancy in the Graph and the total simulation time
    """
    # http://www.regular-expressions.info/floatingpoint.html
    reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')

    lastlines = s.check_output(['tail', '-2', tfile]).strip().split("\n")
    if not reg_flt.match(lastlines[0]):
      raise ValueError('Cannot parse simulation output', tfile)
    else :
      time = float(lastlines[0].split()[0])
      iterations = int(lastlines[-1].split()[-1])
      tot_occ =sum(map(float, lastlines[0].split()[1:]))
      for e, occu in enumerate(lastlines[0].split()[1:]) :
        ss = sorted_nodes[e][0]
        self.node[ss]['occupancy'] = float(occu)/tot_occ

    return time, iterations

  #   prune(self, nodes = None, cutoff = 0.01)
  def prune(self, sorted_nodes, cutoff = None) :
    """ Delete nodes or report them as still reachable. 

    Args:
      sorted_nodes (list): All nodes used in the simulation, sorted by energy

    TODO:
      go through all nodes and update when an inactive node was seen last...
      .. so that we can throw out nodes after some time...
    
    """
    if cutoff is None:
      cutoff = self._cutoff
  
    deleted_nodes = 0
    still_reachables = 0
  
    for ni, data in reversed(sorted_nodes) :
      #print ni, data
      if data['occupancy'] < cutoff :
  
        nbrs = filter(lambda x: self.node[x]['active'], sorted(self.successors(ni), 
            key=lambda x: self.node[x]['energy'], reverse=False))
        best, been = nbrs[0], self.node[nbrs[0]]['energy']
  
        if been > data['energy'] :
          still_reachables += 1
          continue
        
        multibest = {best : been}
      
        (transfer, minbar) = (best, None)
        for e, nbr in enumerate(nbrs[1:]) :
          for mb in multibest.keys():
            always_true = self.add_transition_edge(nbr, mb, ts=ni)
            msE = self.get_saddle(nbr, mb)
            if minbar is None or minbar < (msE - self.node[ni]['energy']):
              (transfer, minbar) = (nbr, msE - self.node[ni]['energy'])
          multibest[nbr] = self.node[nbr]['energy']
          if always_true is False :
            raise TrafoAlgoError('Did not add the transition edge!')
  
        if True: # Set to 'False' to keep all nodes
          self.node[ni]['active']=False
          self.node[ni]['last_seen']=1
          self.node[transfer]['occupancy'] += self.node[ni]['occupancy']
          self.node[ni]['occupancy']=0.0
          deleted_nodes += 1
  
    return deleted_nodes, still_reachables

  def print_to_file(self):
    pass

  def sorted_trajectories_iter(self, sorted_nodes, tfile):
    """ Yields the time course information using a treekin output file.
    
    Args:
      sorted_nodes (list): a list of nodes sorted by their energy
      tfile (str): treekin-output file name.
  
    Yields:
      list: ID, time, occupancy, structure, energy
    """
    # http://www.regular-expressions.info/floatingpoint.html
    reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')
  
    ttime = self._total_time
  
    with open(tfile) as tkn:
      # this is ugly, but used to check if we're at the last line
      prevcourse = []
      tknlines = tkn.readlines()
      for line in tknlines:
        if reg_flt.match(line) :
          course = map(float, line.strip().split())
          time = course[0]
  
          for e, occu in enumerate(course[1:]) :
            # is it above visibility threshold?
            ss = sorted_nodes[e][0]
            sss = ss[0:self._transcript_length]
  
            yield self.node[ss]['identity'], ttime+time, occu, \
                sss, self.node[ss]['energy']
          prevcourse = course
    return 

  def plot(self):
    pass

  def sorted_nodes(self):
    # filter active nodes
    # return sorted list
    sorted_nodes = filter(lambda (n,d): d['active'], 
        sorted(self.nodes(data=True), 
          key=lambda x: x[1]['energy'], reverse=False))


def open_breathing_helices(seq, ss, free = 6):
  """ open all breathable helices, i.e. those that share a base-pair
    with an exterior loop region 
  """
  nbrs = set()
  pt = ril.make_pair_table(ss, base=0)

  # mutable secondary structure 
  nbr = list(ss)

  rec_fill_nbrs(nbrs, ss, nbr, pt, (0, len(ss)), free)

  nbrs.add(''.join(nbr))

  return nbrs

def rec_fill_nbrs(nbrs, ss, mb, pt, (n, m), free):
  """ recursive helix opening
  TODO: Test function, but looks good

  :param nbrs: a set of all neighboring conformations
  :param ss: reference secondary structure
  :param mb: a mutable version of ss, which, after the final round will have
    all breathing helices opened
  :param pt: pair table (zero based)
  :param (n,m): the range of the pt under current investigation
  :param free: number of bases that should be freed

  :return: 
  """
  skip = 0 # fast forward in case we have deleted stuff
  for i in range(n, m) :
    j = pt[i]
    if j == -1 : continue
    if i < skip: continue

    nb = list(ss)
    [o,l] = [0,0]
    [p,q] = [i,j]

    add = True
    while p < q and (l == 0 or o < free):
      if pt[p] != q or p != pt[q] :
        """ this is a multiloop """
        # i,j = 1, len(pt)
        rec_fill_nbrs(nbrs, ''.join(nb), mb, pt, (p,q), free-o)
        add = False
        break

      # remove the base-pairs
      pt[p] = pt[q] = -1
      nb[p] = nb[q] = '.'
      mb[p] = mb[q] = '.'
      o += 2 # one base-pair deleted, two bases freed

      l = 0 # reset interior-loop size
      while (p < q and pt[p+1] == -1):
        [p, l] = [p+1, l+1]
      p += 1
      while (p < q and pt[q-1] == -1):
        [q, l] = [q-1, l+1]
      q -= 1
      o += l

    if add :
      nbrs.add(''.join(nb))
    skip = j+1

  return 

def fold_exterior_loop(md, seq, con, ext_moves):
  """ Constrained folding of the exterior loop.

  All constrained helices are replaced with the motif:
    NNNNNNN
    ((xxx))
  for example a helix with the closing-stack CG-UG:
    CG ~ UG -> CGNNNUG
    (( ~ )) -> ((xxx))
  This reduces the sequence length (n) and therefore the runtime O(n^3),
  and it enables the identification of independent structures with the same
  exterior loop features.
  
  Args:
    md (RNA.md()):      ViennaRNA model details (temperature, noLP, etc.)
    seq (str):          RNA sequence
    con (str):          RNA structure constraint
    ext_moves (dict()): Dictionary storing all mappings from exterior-loop
                        constraints (features) to parents.

  Returns:
    (str, str): 
  """

  spacer = 'NNN'
  pt = ril.make_pair_table(con, base=0)
  ext_seq = ''
  ext_con = ''

  # shrink the sequences
  skip = 0
  for i, j in enumerate(pt):
    if i < skip : continue
    if j == -1 : 
      ext_seq += seq[i]
      ext_con += '.'
    else :
      ext_seq += seq[i] + seq[i+1]
      ext_seq += spacer
      ext_seq += seq[j-1] + seq[j]
      ext_con += '(('
      ext_con += 'x' * len(spacer)
      ext_con += '))'
      skip = j+1

  # If we have seen this exterior loop before, then we don't need to 
  # calculate again, and we have to trace back if the parents are connected.
  if ext_seq in ext_moves :
    css = ext_moves[ext_seq][1]
  else :
    fc_tmp = RNA.fold_compound(ext_seq, md)
    fc_tmp.constraints_add(ext_con, RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP)
    css, cfe = fc_tmp.mfe()
    ext_moves[ext_seq] = [set(), css]
    del fc_tmp
  
  # replace characters in constraint
  c, skip = 0, 0
  for i, j in enumerate(pt):
    if i < skip : continue
    if j == -1 : 
      con = con[:i] + css[c] + con[i+1:]
      c += 1
    else :
      c += len(spacer) + 4
      skip = j+1
  ss = con
    
  return ss, ext_seq

