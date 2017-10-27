
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
  def __init__(self, fullseq, vrna_md):
    """ Conformation Graph initialization.
    """
    super(ConformationGraph, self).__init__()

    self._full_sequence = fullseq
    self._model_details = vrna_md
    self._fold_compound = RNA.fold_compound(fullseq, vrna_md)

    # Adjust simulation parameters
    self._RT = 0.61632077549999997
    if vrna_md.temperature != 37.0 :
      kelvin = 273.15 + vrna_md.temperature
      self._RT = (self._RT/310.15)*kelvin

    # Private instance variables:
    self._transcript_length = 0
    self._total_time = 0
    self._nodeid = 0

    # Default parameters:
    self._cutoff = 0.01
    self._fpath = 20
    self._maxdG = 0 # switched off
    self._k0 = 2e5

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

  def graph_copy(self):
    copy = ConformationGraph(self.full_sequence, self._model_details)
    copy.add_nodes_from(self.nodes(data=True))
    copy.add_edges_from(self.edges(data=True))
    return copy

  def sorted_nodes(self, descending = False):
    """ Returns active nodes sorted by energy. 

    Args:
      data (bool, optional): networkx parameter: return a dictionary of attributes
      descending (bool, optional): sorting parameter. 
        True: energetically high to low. 
        False: energetically low to high.
        Defaults to False.
    
    """
    active = filter(lambda (n,d): d['active'], self.nodes(data=True))
    return sorted(active, key=lambda x: (x[1]['energy'], x[0]), reverse=descending)

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

  def add_transition_edge(self, s1, s2, ts = None, fpath = None, maxdG = None):
    """ Calculates transition rates from direct path barrier heights.

    Uses the *findpath* direct path heuristic to find the lowest energy barrier
    between two secondary structures. Typically s2 is the new, energetically
    better structure, but this is not enforced. In case a transient structure
    "ts" is specified, the rate is computed as the minimum between the direct
    path s1 -> s2 and the indirect folding path s1 -> ts -> s2.
    Rates are computed using the Metropolis model. k0 and RT are global 
    variables in the ConformaitonGraph object.

    Args:
      s1 (str): pathway start secondary structure (must be part of the graph already)
      s2 (str): pathway end secondary structure (may be added to the graph)
      ts (str, optional): transient secondary structure that is not on the
        direct folding path. If a transient structure is specified, the direct
        path barriers s1 -> ts and ts -> s2 must be known already.
      fpath (int, optional): the search width parameter for the findpath routine.
        Defaults to None: using global ConformationGraph parameter.
      maxdG (flt, optional): maximum barrier height. Serves as a cutoff, if the energy
        barrier is greater than maxdG, the edge is not added. 
        maxdG = 0 switches this parameter off.
        Defaults to None: using global ConformationGraph parameter.

    Returns:
      bool: True if the transition edge was added to the graph, False otherwise.
    """

    fullseq = self._full_sequence
    md = self._model_details
    fc = self._fold_compound

    _RT = self._RT
    k0 = self._k0

    if maxdG is None:
      maxdG = self._maxdG

    if fpath is None:
      fpath = self._fpath

    assert s1 != s2
  
    # Lookup the in-direct path barrier first
    if ts : # then we know that the indirect path has to be part of saddles
      tsE1 = self.get_saddle(s1,ts)
      tsE2 = self.get_saddle(s2,ts)
      tsE = max(tsE1, tsE2)

    saddleE = self.get_saddle(s1,s2)

    # Now this is the computationally heavy part ...
    if saddleE is None :
      #if ts:
      #  # Unfortunately, find_path_once returns MAX_INT if no path has been found.
      #  sE = float(fc.find_path_once(s1, s2, int(round(tsE * 100)), fpath))/100
      #  saddleE = min(tsE, sE)
      #  sE = float(fc.find_path_once(s2, s1, int(round(saddleE * 100)), fpath))/100
      #  saddleE = min(saddleE, sE)
      #else :
        saddleE = float(fc.path_findpath_saddle(s1, s2, fpath))/100

    if ts :
      saddleE = min(saddleE, tsE)
 
    # invalid neighbor saddle energies are not stored for future calls.
    if maxdG :
      valid = (ts is not None or saddleE-self.node[s1]['energy'] <= maxdG)
    else :
      valid = True
  
    if valid : # Add the edge.
      e1 = self.node[s1]['energy']
      e2 = self.node[s2]['energy'] if self.has_node(s2) else round(fc.eval_structure(s2), 2)
  
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
    """ Find new secondary structures and add them to the ConformationGraph

    The function supports two move-sets: 1) The mfe structure for the current
    sequence length is connected to all present structures, 2) The conformation
    graph is expanded using helix-breathing. 

    Args:
      extend (int, optional): number of nucleotide extensions before graph
        expansion (updates the global variable transcript length). Defaults to 1.
      exp_mode (str, optional): choose from "mfe-only": only use current mfe
        structure as potential new neighbor. "breathing-only": only use breathing
        neighborhood. "default": do both mfe and breathing.
      mfree (int, optional): minimum number of freed bases during a
        helix-opening step. Defaults to 6.
      cutoff (flt, optional): Occupancy cutoff for neighbor generation.
        Defaults to None: using global ConformationGraph parameter.

    Returns:
      int: Number of new nodes
    """
    if cutoff is None:
      cutoff = self.occupancy_cutoff

    fseq = self.full_sequence
    self._transcript_length += extend
    if self._transcript_length > len(fseq):
      self._transcript_length = len(fseq)
    seq = self.transcript

    csid = self._nodeid
    md = self._model_details
    fc_full = self._fold_compound

    if exp_mode not in ['default', 'mfe-only', 'breathing-only']:
      raise TrafoUsageError('unknown expansion mode')

    # Calculate MFE of current transcript
    fc_tmp = RNA.fold_compound(seq, md)
    ss, mfe = fc_tmp.mfe()
    future = '.' * (len(fseq)-len(seq))
    ss = ss + future
    #print "{}\n{} {:6.2f}".format(seq, ss, mfe)

    # If there is no node because we are in the beginning, add the node,
    # otherwise, try to add transition edges from every node to MFE.
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

      # Initialize a dictionary to store the feature expansion during each
      # graph expansion round: ext_moves[ext_seq] = [set((con,paren),...),
      # structure] where ext_seq = exterior-loop sequence with ((xxx))
      # replacing constrained elements
      ext_moves = dict()
      # Every neighbor generation can only produce energetically better
      # neighbors, so they are sorted upfront. Note that there are inactive
      # nodes that might become active during expansion, but their occupancy
      # will not change.
      for ni, data in sorted(self.nodes(data=True), 
          key=lambda x: x[1]['energy'], reverse=True):
        if data['active'] == False : continue
        en  = data['energy']
        occ = data['occupancy']
        if occ < cutoff : continue

        # short secondary structure (without its future)
        sss = ni[0:len(seq)]

        # compute a set of all helix breathing open steps
        opened = open_breathing_helices(seq, sss, mfree)

        # do a constrained exterior loop folding for all of them and then
        # connect them to the present conformation "ni", but also to each
        # other.
        connect = [ni]
        for onbr in opened :
          nbr, ext_seq = fold_exterior_loop(md, seq, onbr, ext_moves)

          future = '.' * (len(ni) - len(nbr))
          nbr += future

          if ni == nbr: 
            continue

          if self.has_edge(ni, nbr):
            self.node[nbr]['active'] = True # in case it was there but inactive
            continue

          for con in connect:
            if con == nbr: 
              continue
            if self.has_node(nbr):
              if self.add_transition_edge(con, nbr):
                self.node[nbr]['active'] = True # in case it was there but inactive
            elif self.add_transition_edge(con, nbr) :
              enbr = round(fc_full.eval_structure(nbr), 2)
              self.node[nbr]['energy'] = enbr
              self.node[nbr]['active'] = True
              self.node[nbr]['last_seen'] = 0
              self.node[nbr]['occupancy'] = 0.0
              self.node[nbr]['identity'] = self._nodeid
              self._nodeid += 1
            else :
              # Replace this with a continue, or logging... but first see some examples.
              # raise DebuggingAlert("# helix breathing: could not add transition edge!")
              continue

          if self.has_node(nbr) and self.node[nbr]['active'] :
            connect.append(nbr)

          # Now connect the neighbor with *historic* transitions of parents.
          # We store the exterior-open neighbor here, that means there are
          # three possible reasons for duplication:
          #   1) different (or longer) helix was opened / same historic features
          #   2) the same helix was opened / difference is in historic features
          #   3) different helix / different history
          if ext_moves[ext_seq][0] :
            for (parent, child) in ext_moves[ext_seq][0] :
              assert parent != ni # Parents may never be the same
              if child == nbr : 
                # the parents differed in breathing helices, no historic differences
                continue 

              if self.has_edge(parent, ni) :
                if self.has_node(child) and self.has_node(nbr):
                  if self.has_edge(nbr, child): 
                    continue

                  # TODO: Calculate saddleE from saddleE of parents?
                  sP = self.get_saddle(ni, parent)
                  sC1 = round(self.node[child]['energy'] + sP - self.node[parent]['energy'], 2)
                  sC2 = round(self.node[nbr]['energy'] + sP - self.node[ni]['energy'], 2)

                  if sC1 == sC2: # avoid findpath!
                    self.add_weighted_edges_from([(nbr, child, None)])
                    self.add_weighted_edges_from([(child, nbr, None)])
                    self[nbr][child]['saddle'] = sC1
                    self[child][nbr]['saddle'] = sC1
                    if self.add_transition_edge(nbr, child, maxdG = 0):
                      self.node[nbr]['active'] = True # in case it was there but inactive
                      self.node[child]['active'] = True # in case it was there but inactive
                    else :
                      self.remove_edge(nbr,child)
                      self.remove_edge(child,nbr)
                      # Replace this with a continue, or logging... but first see some examples.
                      raise DebuggingAlert("# helix breathing: did not add historic edge! It is safe to replace this with a continue, but tell stefan about it!")
                  else:
                    if self.add_transition_edge(nbr, child):
                      self.node[nbr]['active'] = True # in case it was there but inactive
                      self.node[child]['active'] = True # in case it was there but inactive
                else :
                  # TODO: Need to think about this more once a few examples come up ...
                  raise DebuggingAlert("# helix breathing: should we add the edge?! It is safe to replace this with a continue, but tell stefan about it!")
          # Track the final structure, every new identical ext-change will be
          # connected, if the parents were connected.
          ext_moves[ext_seq][0].add((ni, nbr))

    if not self.has_node(ss) or (not self.node[ss]['active']):
      print "# WARNING: ", ss, "[mfe secondary structure not connected]"

    # Post processing of graph after expansion: 
    #   remove nodes that have been inactive for a long time.
    for ni in self.nodes() :
      if self.node[ni]['active'] == False :
        self.node[ni]['last_seen'] += 1
      else :
        self.node[ni]['last_seen'] = 0
      if self.node[ni]['last_seen'] >= 5:
        self.remove_node(ni)

    return self._nodeid-csid

  def coarse_grain(self, minh = 1):
    """ Graph condensation step. 

    Processes an energetically sorted list of structures (high to low energies)
    and tries to merge their occupancy into a neighboring, better conformation
    (thus forming a macro-state). The current structure can be merged (and
    therefore removed from the graph) if there exists a neighbor that has
    better energy and the energy barrier is lower than the *minh* parameter. If
    there are multiple better neighbors with the same transition energy
    barrier, then the occupancy is transfered to the neighbor with the lowest
    energy, in the degenerate case the lexicographically first structure is
    chosen.

    Args:
      minh (flt, optional): Minimum energy barrier between structures.

    Returns:
      dict[node] = node: A mapping from deleted nodes to macro-state
    """

    merged_nodes = dict()
    merged_to = dict()

    # sort by energy (high to low)
    for ni, data in sorted(self.nodes(data=True), key=lambda x: (x[1]['energy'], x), reverse=True):
      if data['active'] == False : continue
      en  = data['energy']

      # get all active neighbors (low to high)
      nbrs = filter(lambda x: self.node[x]['active'], 
          sorted(self.successors(ni), key=lambda x: (self.node[x]['energy'], x), reverse=False))

      if nbrs == []:
        break

      # lowest neighbor structure and energy
      best, been = nbrs[0], self.node[nbrs[0]]['energy']
  
      if been - en > 0.0001:
        # still reachable, don't remove this node ... 
        continue

      # among all energetically better neighbors, find the neighbor with the
      # lowest energy barrier ...
      (transfer, minsE) = (best, self.get_saddle(ni, best))
      for e, nbr in enumerate(nbrs[1:]) :
        if self.node[nbr]['energy'] - en >= 0.0001 :
          break
        sE = self.get_saddle(ni, nbr)
        if sE - minsE < 0.0001 :
          (transfer, minsE) = (nbr, sE)
      
      if minsE - en - minh > 0.0001 : # avoid precision errors
        # do not merge, if the barrier is too high.
        continue

      # connect all neighboring nodes with each other
      for e, nb1 in enumerate(nbrs, 1) :
        for nb2 in nbrs[e:] :
          always_true = self.add_transition_edge(nb2, nb1, ts=ni)
          if always_true is False :
            raise TrafoAlgoError('Did not add the transition edge!')
  
      # remove the node
      self.node[ni]['active']=False
      self.node[ni]['last_seen']=1
      self.node[transfer]['occupancy'] += self.node[ni]['occupancy']
      self.node[ni]['occupancy']=0.0

      merged_nodes[ni] = transfer
      if transfer in merged_to :
        merged_to[transfer].append(ni)
      else :
        merged_to[transfer] = [ni]

      if ni in merged_to :
        fathers = merged_to[ni]
        for f in fathers:
          merged_nodes[f] = transfer
          merged_to[transfer].append(f)
        del merged_to[ni]

    return merged_nodes

  def simulate(self, t0, t8, tmpfile=None): 
    # treekin wrapper function using:
    #   "self.get_simulation_files_tkn" 
    #   "self.update_occupancies_tkn"
    raise NotImplementedError

  def get_simulation_files_tkn(self, name, binrates = True):
    """ Print a rate matrix and the initial occupancy vector.

    This function prints files and parameters to simulate dynamics using the
    commandline tool treekin. A *.bar file contains a sorted list of present
    structures, their energy and their neighborhood and the corresponding
    energy barriers. A *.rts or *.rts.bin file contains the matrix of
    transition rates either in text or binary format. Additionaly, it returns 
    a vector "p0", which contains the present occupancy of structures. The 
    order or elements in p0 contains

    Note: 
      A *.bar file contains the energy barriers to transition between local
      minima. In contrast to files produced by `barriers`, where local minimum
      is always *directly* connected to an energetically better local minimum, 
      here a path towards the MFE structure can proceed via an energetically
      worse structure first. 

    Args:
      name (str): Name of output files name.bar, name.rts, name.rts.bin.  
      binrates (bool, optional): Print rates in binary format or text format.
        Defaults to True: binary format.

    """
    seq = self.transcript

    sorted_nodes = self.sorted_nodes(descending = False)

    bfile = name+'.bar'
    rfile = name+'.rts'
    brfile = rfile+'.bin'
    p0 = []
    
    with open(bfile, 'w') as bar, open(rfile, 'w') as rts, open(brfile, 'w') as brts :
      bar.write("     {}\n".format(seq))
      brts.write(pack("i", len(sorted_nodes)))
      for e, (ni, data) in enumerate(sorted_nodes, 1) :
        # Calculate barrier heights to all other basins.
        nMsE = set()
        for ee, (be, _) in enumerate(sorted_nodes, 1) :
          if e == ee :
            continue
          sE = self.get_saddle(be,ni)
          if sE is not None:
            nMsE.add((ee, sE))
        mystr = ' '.join(map(lambda(x,y):'({:3d} {:6.2f})'.format(x,y-data['energy']), 
            sorted(list(nMsE), key=lambda x:x[0])))

        # Print structures and neighbors to bfile:
        bar.write("{:4d} {} {:6.2f} {}\n".format(
          e, ni[:len(seq)], data['energy'], mystr))

        # Add ni occupancy to p0
        if data['occupancy'] > 0 :
          p0.append("{}={}".format(e,data['occupancy']))

        # Print rate matrix to rfile and brfile
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

    return [bfile, brfile if binrates else rfile, p0, sorted_nodes]

  #   update_time_and_occupancies_tkn(self, tfile)
  def update_occupancies_tkn(self, tfile, sorted_nodes):
    """
      Update the occupancy in the Graph and the total simulation time
    """
    # http://www.regular-expressions.info/floatingpoint.html
    reg_flt = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?.')

    lastlines = s.check_output(['tail', '-2', tfile]).strip().split("\n")
    if not reg_flt.match(lastlines[0]):
      raise TrafoAlgoError('Cannot parse simulation output', tfile)
    else :
      time = float(lastlines[0].split()[0])
      iterations = int(lastlines[-1].split()[-1])
      tot_occ =sum(map(float, lastlines[0].split()[1:]))
      for e, occu in enumerate(lastlines[0].split()[1:]) :
        ss = sorted_nodes[e][0]
        self.node[ss]['occupancy'] = float(occu)/tot_occ

    return time, iterations

  def prune(self, cutoff = None, maxh = None) :
    """ Delete nodes or report them as still reachable. 

    Use the occupancy cutoff to choose which nodes to keep and which ones to
    remove. Every node with occuancy < cutoff will be removed and its neighbors
    connected with each other. You may set the *maxh* parameter to reject the 
    removal of a node that has a very high energy barrier to all its neighbors.

    Args:
      cutoff (flt, optional): Occupancy cutoff for neighbor generation.
        Defaults to None: using global ConformationGraph parameter.
      maxh (flt, optional): Don't remove structures that are separated with 
        an energy barrier higher than maxh.

    Returns:
      int, int, int: 
        number of deleted nodes, 
        number of still reachable nodes, 
        number of rejected deletions due to maxh
    
    """
    if cutoff is None:
      cutoff = self._cutoff
  
    deleted_nodes = 0
    still_reachables = 0
    rejected = 0

    for ni, data in self.sorted_nodes(descending=False) : # sort high to low..
      if data['occupancy'] - cutoff > 0.0000001 :
        continue
      en  = data['energy']
  
      # get all active neighbors (low to high)
      nbrs = filter(lambda x: self.node[x]['active'], sorted(self.successors(ni), 
          key=lambda x: self.node[x]['energy'], reverse=False))

      # lowest neighbor structure and energy
      best, been = nbrs[0], self.node[nbrs[0]]['energy']
  
      if been - en > 0.0001:
        still_reachables += 1
        continue
      
      # among *all* neighbors, find the neighbor with the lowest energy barrier
      (transfer, minsE) = (best, self.get_saddle(ni, best))
      for e, nbr in enumerate(nbrs[1:]) :
        sE = self.get_saddle(ni, nbr)
        if sE - minsE < 0.0001 :
          (transfer, minsE) = (nbr, sE)
  
      if maxh and (minsE - en - maxh > 0.0001) : 
        # do not merge, if the barrier is too high.
        rejected += 1
        continue

      # connect all neighboring nodes with each other
      for e, nb1 in enumerate(nbrs, 1) :
        for nb2 in nbrs[e:] :
          always_true = self.add_transition_edge(nb2, nb1, ts=ni)
          if always_true is False :
            raise TrafoAlgoError('Did not add the transition edge!')

      # remove the node
      self.node[ni]['active']=False
      self.node[ni]['last_seen']=1
      self.node[transfer]['occupancy'] += self.node[ni]['occupancy']
      self.node[ni]['occupancy']=0.0
      deleted_nodes += 1
  
    return deleted_nodes, still_reachables, rejected

  def sorted_trajectories_iter(self, sorted_nodes, tfile, softmap = None):
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
  
          # softmap hack:
          # preprocess the timeline by merging all states using the softmap
          if softmap:
            macrostates = [0] * len(course)
            macromap = dict()
            for e, occu in enumerate(course[1:]) :
              ss = sorted_nodes[e][0]
              #print ss, occu

              if ss in softmap : # map occupancy to (energetically better) softmap[ss]
                mapss = softmap[ss]
                mapid = macromap[mapss]
              else : 
                # we *must* have seen this state before, given there are no degenerate
                # sorting errors...
                mapid = e+1
                macromap[ss] = mapid

              macrostates[mapid] += occu

            course[1:] = macrostates[1:]

          for e, occu in enumerate(course[1:]) :
            # is it above visibility threshold?
            ss = sorted_nodes[e][0]
            sss = ss[0:self._transcript_length]
  
            yield self.node[ss]['identity'], ttime+time, occu, \
                sss, self.node[ss]['energy']
          prevcourse = course
    return 


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

