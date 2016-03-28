#!/usr/bin/env python

""" interkin.py
  Written by Stefan Badelt (stef@tbi.univie.ac.at)

  Department of Theoretical Chemistry, University of Vienna
  http://www.tbi.univie.ac.at

  requries: BarMap.py => ViennaRNA-v2.1.9, barriers-v1.6

  main(): interkin.py < file.fasta > graph.dump 2> loginfo

  TODO: Internal ODE solving instead of 'graph.dump'

  vim-config = set: ts=2 et sw=2 sts=2
"""

import re
import os
import sys
import string
import argparse
import pandas as pd
import networkx as nx
import collections as c
import matplotlib.pyplot as plt
import math
#import cofold_plot as conc

import rnaworld as nal
import RNA

def module_exists(module_name):
  try:
    __import__(module_name)
  except ImportError:
    return False
  else:
    return True

def draw_d3_graph(RG):
  if module_exists('d3py'):
    import d3py
    with d3py.NetworkXFigure(RG, width=1500, height=1500) as p:
      p += d3py.ForceLayout()
      p.css['.node'] = {'fill': 'blue', 'stroke': 'magenta', 'title': 'bla'}
      p.css['.link'] = {'stroke': 'red', 'stoke-width': '3px'}
      p.show()
  else :
    print >> sys.stderr, 'Need to install Module: "d3py"'

  return

def get_species(fullseq, homophobic) :
  ''' Split a sequence of multiple molecules ('&') into a list:
  :param: homophobic Choose whether to return homo-dimers or not

  :return: [(name, seq), (name&name, seq&seq), ...]
  '''
  subseqs  = c.OrderedDict.fromkeys(fullseq.split('&'))
  species  = list(string.uppercase[:len(subseqs)])
  monomers = zip(species, subseqs)

  dim = len(monomers)
  dimers = []
  for x in range(dim):
    for y in range(x, dim):
      if (homophobic and x == y) : continue
      (speX, seqX) = monomers[x]
      (speY, seqY) = monomers[y]
      dimers.append((speX+'&'+speY, seqX+'&'+seqY))
  return monomers + dimers

def dump_graph(graph, nodes, args):
  ''' Dump info of the graph into a file:
    translate that file into Sundials CVODE using:
    cat file | /path/to/sundials/scripts/tryit.pl -n simulation

    TODO: try 'to_edgelist'
  '''
  gfile = args.name+'_graph.dump'
  with open(gfile, 'w') as g :
    for n in nodes:
      if nodes[n] :
        g.write("{:s}\n".format(':'.join(list(n))))
    for n, nbrs in graph.adjacency_iter():
      n = ':'.join(n)
      for nbr, edict in nbrs.items():
        #if len(edict) > 1 : print "Multiedge!!!", n, nbr, edict
        nbr = ':'.join(nbr)
        for d in edict.values() :
          #print("%s %s %.3f\n" % (n,nbr,d['weight']))
          g.write("%s %s %.3f\n" % (n,nbr,d['weight']))

def all_forward_hyper_rates(RG, (nodeA, nodeB)):
  """ forward rate: A + B => AB """
  #print nodeA, "x", nodeB
  forward = []
  for suc1 in RG.successors_iter(nodeA) :
    if suc1[0] != 'HYPER' : continue
    for suc2 in RG.successors_iter(nodeB) :
      if suc1 is not suc2 : continue
      edict = RG.get_edge_data(nodeA, suc1)
      rate = edict[0]['weight']
      forward.append((rate))
  return forward

def all_backward_hyper_rates(RG, (nodeAB,)):
  """ backward rate: AB => A + B """
  #print nodeAB
  backward = []
  for succ in RG.successors_iter(nodeAB) :
    if succ[0] != 'HYPER' : continue
    edict = RG.get_edge_data(nodeAB, succ)
    rate = edict[0]['weight']
    backward.append((rate))
  return backward

def rate_balance(RG, nprobs):
  """ detailed balance = sum over all Pi*kij - Pj*kji """
  detbal = []
  for node in nprobs :
    for succ in RG.successors_iter(node):
      """ Sum over all Pi * kij """
      if succ[0] == 'HYPER' : continue
      edict = RG.get_edge_data(node, succ)
      rate = edict[0]['weight']
      fw = nprobs[node]*rate
      edict = RG.get_edge_data(succ, node)
      rate = edict[0]['weight']
      bw = nprobs[succ]*rate
      if False and fw-bw > 0.1 :
        print node, succ, fw-bw
      detbal.append(abs(fw-bw)/2)
  return sum(detbal)

def add_basins(RG, bfile, spe, seq, nodelist, pure):
  """ Read the basin free energies from barrier tree files 
    Add them up to get the total partition function Z
    store the basin partition functions in the array zs
    check detailed balance and compute reaction constants
  """
  homodimer = ('&' in spe and spe.split('&')[0] == spe.split('&')[1])
  kT=0.61632077549999997
  Z = 0.
  zs = []
  nodes = []
  for line in nal.parse_barfile(bfile, seq=seq) :
    if nodelist and nodelist[(spe, line[1])] == 0 : 
      """ this is only needed in the bionly case """
      rescue_flag = 1
      for cplx in get_complexes((spe, line[1])) :
        if nodelist[cplx] == 0 :
          rescue_flag = 0
          break
      if not rescue_flag :
        """ normally minima high up in the cofold tree, because
        truemin tree has additional minima """
        #print "lost:", spe, line
        continue
    if '&' in spe and pure != 0 :
      """ Filter the cofolded landscape """
      cplx = get_complexes((spe, line[1]))
      if len(cplx) > 1 :
        """ if it is a homodimer: zA*zA*2 """
        if pure == 1 or pure == 3 :
          Gb = float(line[-1])
          z = math.e**(-Gb/kT)
          Z += z
          zs.append(z)
          nodes.append(cplx)
      else :
        """ if it is a homodimer: zAA*2 """
        if pure == 2 or pure == 3 :
          Gb = float(line[-1])
          z = math.e**(-Gb/kT)
          Z += z
          zs.append(z)
          nodes.append(cplx)
    else :
      """ unimolecular stuff, only remember the corresponding
       vertices, do the balance check lateron """
      Gb = float(line[-1])
      z = math.e**(-Gb/kT)
      Z += z
      zs.append(z)
      nodes.append((spe, line[1]))

  """ symmetry correction, not considered in kinetic simulations! """
  if False and homodimer :
    print "correcting"
    Z /= 2
    zs = [x/2 for x in zs]

  """ 
    Translate grandient basins partition functions into probabilities and 
    Z into the ensemble free energy 
  """
  zs = [x/Z for x in zs]
  Z = -kT*math.log(Z)

  """ compute the forward and backward rates for each node in the graph """
  detbal = 0
  if pure == 1 :
    for cplx, prob in zip(nodes, zs) :
      for rate in all_forward_hyper_rates(RG, cplx) :
        detbal += rate * prob
  elif pure == 2 :
    for cplx, prob in zip(nodes, zs) :
      for rate in all_backward_hyper_rates(RG, cplx) :
        detbal += rate * prob
  elif pure == 0 :
    nprobs = c.OrderedDict(zip(nodes,zs))
    detbal = rate_balance(RG, nprobs)
  elif pure == 3 :
    nprobs = c.OrderedDict()
    for nlist, prob in zip(nodes, zs) :
      for n in nlist :
        nprobs[n]=prob
    detbal = rate_balance(RG, nprobs)

  cutoff = 0.01
  for e, p in enumerate(zs) :
    if p < cutoff :
      zs = zs[:e]
      break
  if detbal == 0 :
    detbal = 1
  return Z, zs, detbal

'''
def basin_free_energies(
  tmpcut = RNA.cvar.cut_point

  [seqA, seqB] = seq.split('&')

  RNA.cvar.cut_point = len(seqA)+1
  [stmp,G_A,G_B,G_AB,waste]   = RNA.co_pf_fold(seqA+seqB,None)

  if (seqA == seqB) :
    G_AA = G_BB = G_AB
  else :
    [stmp,tmp1,tmp2,G_AA,waste] = RNA.co_pf_fold(seqA+seqA,None)
    RNA.cvar.cut_point = len(seqB)+1
    [stmp,tmp1,tmp2,G_BB,waste] = RNA.co_pf_fold(seqB+seqB,None)

  RNA.cvar.cut_point = tmpcut
'''

def bimolecular_cofold_basins(RG, spe, seq, args, ns, bionly):
  ''' Calculate ensemble free energies similar to co_pf_fold() 
    Do it for 4 landscapes:
    G_A : Monomers of type A
    G_B : Monomers of type B
    G_AB: Dimers AB *patch*
    G_AB_TOT: Monomers of A, B and Dimers AB
  '''
  name = args.name
  args.name = name+"_"+spe
  args.name = args.name.translate(string.maketrans("&", "_"))

  [bfile, efile, rfile, ps] = bm.barriers(seq, args)

  #[sfile, bfile, efile, rfile, psfile] = nal.sys_barriers(args.name, seq, args.name+'.spt', 
  #      barriers=args.barriers,
  #      minh=args.b_minh,
  #      maxn=args.b_maxn,
  #      k0=args.k0,
  #      temp=args.temperature,
  #      noLP=args.noLP,
  #      moves='single-base-pair',
  #      gzip=True,
  #      rates=True,
  #      bsize=False,
  #      saddle=False,
  #      mfile=mfile,
  #      force=args.force,
  #      verb=args.verbose)

  """ get the fake-dimer basins, homo-double counts removed """
  Gm, Mprobs, fw = add_basins(RG, bfile, spe, seq, ns, 1)

  """ get the true-dimer basins, homo-symmetry corrected """
  Gd, Dprobs, bw = add_basins(RG, bfile, spe, seq, ns, 2)
    
  """ get Z_tot = Z_A * Z_B + Z_AB """
  if bionly :
    Gt, Tprobs, detbal = add_basins(RG, bfile, spe, seq, ns, 3)
    print "Cofold detailed blance != 0 ? ", detbal

  if False :
    """ Now lets see which species should be populated """
    print Mprobs
    print Dprobs
    #print Tprobs

  print "K from bimolecular rates:", fw, bw, "=", fw/bw
  args.name = name

  return [Gm,Gd]

def unimolecular_basins(RG, ispe, iseq, args, ns):
  ''' Calculate ensemble free energies similar to co_pf_fold() 
    Do it for 4 landscapes:
    G_A : Monomers of type A
    G_B : Monomers of type B
    G_AB: Dimers AB *patch*
    G_AB_TOT: Monomers of A, B and Dimers AB
  '''

  species = zip(ispe.split('&'), iseq.split('&'))
  species.append((ispe,iseq))

  G = [] # A, B, AB
  for (spe, seq) in species :
    name = args.name
    args.name = name+"_"+spe
    if '&' in args.name :
      """ check true dimers """
      args.name = args.name.translate(string.maketrans("&", "_"))
      args.spatch.append('-d')
      [bfile, efile, rfile, ps] = bm.barriers(seq, args)
      Gd, probs, detbal = add_basins(RG, bfile, spe, seq, ns, 0)
      print "Detailed balance != 0 ? ", detbal
      args.spatch.pop()
      G.append((Gd))
    else :
      [bfile, efile, rfile, ps] = bm.barriers(seq, args)
      Gm, probs, detbal = add_basins(RG, bfile, spe, seq, ns, 0)
      print "Detailed balance != 0 ? ", detbal
      G.append((Gm))
    args.name = name
  return G

def interkin_equilibrium(species, args):
  ''' Check if the gradient basins from different barrier-trees match
    * assess how big the error from coarse-graining L_AB is
    * compute the new equilibrium distribution
  '''
  verb = args.verbose
  kT=0.61632077549999997
  TreeData = c.defaultdict(list)
  monomers = []
  fakedimers = []
  truedimers = []
  # Get all the barrier-trees and store the info in TreeData
  for (spe, seq) in species :
    patchdata = []
    name = args.name
    args.name = name+"_"+spe
    if '&' in args.name :
      args.name = args.name.translate(string.maketrans("&", "_"))

      ''' Do an additional run first: 'true-dimers-only' '''
      args.spatch.append('-d')
      [bfile, efile, rfile, ps] = bm.barriers(seq, args)
      patchdata  = nal.parse_barfile(bfile, seq=seq)
      args.spatch.pop()
    [bfile, efile, rfile, ps] = bm.barriers(seq, args)
    bardata = nal.parse_barfile(bfile, seq=seq)
    args.name = name

    ''' Write the Barriers Output into TreeData dictionary '''
    for line in bardata :
      TreeData[(spe, line[1])] = line
      if '&' in spe :
        cplx = get_complexes((spe, line[1]))
        if len(cplx) > 1 :
          fakedimers.append(cplx)
        else :
          truedimers.append(cplx)
      else :
        monomers.append((spe, line[1]))
    for line in patchdata :
      cplx = get_complexes((spe, line[1]))
      if len(cplx) > 1 :
        sys.exit('something went terribly wrong')
      else :
        newid = list(cplx[0])
        newid[0] = '&'+newid[0]+'&'
        TreeData[tuple(newid)] = line

  ''' Compare coarse-grainings of uni- and bi-molecular landscapes 
    In the united landscape L_{AB}, gradient basins can be messed up:
    * Gradient walks starting in a dimer conformation can end in a basin that
      contains two monomers. Hence, the monomer-basins can be bigger than they
      actually should be!
    * At the high parts of the enery landscapes, monomer-basins can also be 
      smaller, as the number of suboptimals limits the structures contributing
      to a gradient basin. 
  '''
  ''' Z_AB = Z_A * Z_B <==> G_AB = G_A + G_B '''
  spicy, fakeG, trueG = [], [], []
  for cplx in fakedimers : 
    """compare two monomers with fakedimer"""
    bimol = [x+'&'+y for (x,y) in zip(*cplx)]
    monoPF = []
    for cx in cplx:
      monoPF.append(float(TreeData[cx][-1]))
    #print tuple(bimol), "fake:", float(TreeData[tuple(bimol)][-1]), \
    #    "true:", sum(monoPF), monoPF
    spicy.append(tuple(bimol))
    """ original tree data from cofold """ 
    fakeG.append(float(TreeData[tuple(bimol)][-1]))
    trueG.append(sum(monoPF))

  diffG = [x-y for x,y in zip(fakeG, trueG)]
  print "Total difference in bfe beween true and cofold-dimers:", sum(diffG)
  F_ApB = -kT*math.log(sum([math.e**(-x/kT) for x in fakeG]))
  T_ApB = -kT*math.log(sum([math.e**(-x/kT) for x in trueG]))
  D_ApB = F_ApB + sum(diffG) #-kT*math.log(sum([math.e**(-x/kT) for x in diffG]))
  print "F", F_ApB 
  print "T", T_ApB 
  print "D", D_ApB 

  ''' Z'_AB = Z_AB <==> G'_AB = G_AB '''
  spicy, fakeG, trueG = [], [], []
  for cplx in truedimers :
    'compare truedimers (patch) with fakedimer'
    bimol = cplx[0]
    newid = list(cplx[0])
    newid[0] = '&'+newid[0]+'&'
    if TreeData[tuple(newid)] :
      #print bimol, "fake:", TreeData[tuple(bimol)][-1], \
      #    "true:", TreeData[tuple(newid)][-1] 
      spicy.append(tuple(bimol))
      """ original tree data from cofold """ 
      fakeG.append(float(TreeData[tuple(bimol)][-1]))
      """ tree data from patch """ 
      trueG.append(float(TreeData[tuple(newid)][-1]))
    else :
      """ Not found in True landscape """
      Gb = float(TreeData[tuple(bimol)][-1])

  #print [x-y for x,y in zip(fakeG, trueG)]
  diffG = [x-y for x,y in zip(fakeG, trueG)]
  print "Total difference in bfe beween true and cofold-dimers:", sum(diffG)
  F_AB = -kT*math.log(sum([math.e**(-x/kT) for x in fakeG]))
  T_AB = -kT*math.log(sum([math.e**(-x/kT) for x in trueG]))
  D_AB = F_AB + sum(diffG) #-kT*math.log(sum([math.e**(-x/kT) for x in diffG]))
  print "F", F_AB 
  print "T", T_AB 
  print "D", D_AB 

  a0s = map(lambda x: 10**-x, range(9,8,-3))
  b0s = map(lambda x: 10**-x, range(6,1,-1))

  K = math.e**((D_ApB-D_AB)/kT)
  eqdata = conc.bimol_equilibrium(K, a0s, b0s)
  [a0, b0, eA, eB, eAB] = zip(*eqdata)
  conc.bimol_plot(eqdata, 'diff_test_plot')

  for cplx in monomers :
    pass

  return


""" Folding Kinetics using the rnaworld library """

def cofold_barriers(_name, species,
    s_opts = ['|', 'pylands_spatch.py'],
    s_ener = None,
    s_maxn = 59000000,
    barriers='barriers',
    minh=0.001,
    maxn=50,
    k0=1.0,
    temp=37.0,
    moves='single-base-pair',
    gzip=True,
    rates=True,
    bsize=False,
    saddle=False,
    bionly=False,
    force=False,
    verb=False):
  """ Simulate inter-molecular folding
  Compute folding kinetics 
    returns the biggest connected component and a list of nodes
    TODO: 
      *) @homo-dimers: RNAsubopt-barriers does not correct for two-fold symmetry!
  """
  
  noLP=False # MUST BE FALSE!
  circ=False # MUST BE FALSE!
  mfile=''

  if minh > 0.5 :
    print "Warning: increasing the minh option for cofolded barriers may " + \
    "destroy the separation of monomer and dimer gradient basins!"

  RG = nx.MultiDiGraph()
  sortednodes = c.OrderedDict()
  hypercount  = [0]

  for spe, seq in species :
    name = _name + "_" + spe

    if '&' in name :
      """ true dimer only run """
      name = name.translate(string.maketrans("&", "_"))
      
      if s_ener is None :
        myener, nos = nal.sys_subopt_range(seq, nos=s_maxn, maxe=20, verb=verb)
      else :
        myener = s_ener

      if bionly is False :
        _fname = name + "_truedimers"
        s_opts.append('-d')
        sfile = nal.sys_suboptimals(_fname, seq, 
            ener=myener, opts=s_opts, verb=verb, force=force)
        [sfile, bfile, efile, rfile, psfile] = nal.sys_barriers(_fname, seq, sfile, 
            barriers=barriers, minh=minh, maxn=maxn, k0=k0, temp=temp, 
            verb=verb, force=force)
        # this function modifies sortednodes and hypercount... not nice!
        add_edges(RG, bfile, rfile, 'uni', spe, seq, sortednodes, hypercount)
        s_opts.pop()

      _fname = name
      sfile = nal.sys_suboptimals(_fname, seq, 
          ener=myener, opts=[], verb=verb, force=force)
      [sfile, bfile, efile, rfile, psfile] = nal.sys_barriers(_fname, seq, sfile, 
          barriers=barriers, minh=minh, maxn=maxn, k0=k0, temp=temp, 
          verb=verb, force=force)
      if bionly :
        add_edges(RG, bfile, rfile, 'bionly', spe, seq, sortednodes, hypercount)
      else :
        add_edges(RG, bfile, rfile, 'bi', spe, seq, sortednodes, hypercount)
    elif bionly :
      continue
    else :
      if s_ener is None :
        myener, nos = nal.sys_subopt_range(seq, nos=s_maxn, maxe=20, verb=verb)
      else :
        myener = s_ener
      _fname = name
      sfile = nal.sys_suboptimals(_fname, seq, 
          ener=myener, opts=[], verb=verb, force=force)
      [sfile, bfile, efile, rfile, psfile] = nal.sys_barriers(_fname, seq, sfile, 
          barriers=barriers, minh=minh, maxn=maxn, k0=k0, temp=temp, 
          verb=verb, force=force)
      add_edges(RG, bfile, rfile, 'uni', spe, seq, sortednodes, hypercount)

    print "connected components", 
    print [len(comp) for comp in 
        sorted(nx.strongly_connected_components(RG), key=len, reverse=True)]

  for e, comp in enumerate(
      sorted(nx.strongly_connected_components(RG), key=len, reverse=True)) :
    #print comp
    if e == 0 : continue
    print "Discarding connected component!", len(comp), comp
    for i in comp: sortednodes[i] = 0

  return max(nx.strongly_connected_component_subgraphs(RG), key=len), sortednodes

def add_edges(RG, bfile, rfile, mode, spe, seq, snodes, hcount):
  """ Add the transition edges, no symmetry correction! 
  WARNING: this function modifies snodes and hcount, that is a bit confusing...
  
  """
  RM = []
  lmins = [] #map(lambda(x): (spe, x[1]), nal.parse_barfile(bfile, seq=seq))
  for line in nal.parse_barfile(bfile, seq=seq):
    lmins.append((spe, line[1]))
    if (spe, line[1]) not in snodes :
      snodes[(spe,line[1])] = 0
  with open(rfile) as rates :
    for line in rates :
      """ Can we correct rates for two-fold rotational symmetry ? """
      RM.append((map(float, line.strip().split())))
  pDF = pd.DataFrame(RM, columns=lmins, index=lmins)

  if mode == 'uni' :
    for col in pDF : 
      for row in pDF :
        if row == col : continue
        if pDF[col][row] == 0. : continue
        if pDF[row][col] == 0. :
          print "backward rate missing!!!"
        #The forward rate row => col is in RM[col][row]
        #print "Adding:", row, col, pDF[col][row]
        RG.add_weighted_edges_from([(row, col, pDF[col][row])])
        snodes[row]=1
        snodes[col]=1
  elif mode == 'bi' or mode == 'bionly' :
    #fakemin
    for col in pDF : 
      for row in pDF :
        if row == col : continue
        if pDF[col][row] == 0. : continue
        #The forward rate row => col is in RM[col][row]
        cplR = get_complexes(row)
        cplC = get_complexes(col)
        # good enough to identify pairwise interactions:
        if len(cplR) == len(cplC) : 
          """ WARNING: adding unimolecular rate from cofold tree """
          if mode == 'bi' : continue
          if len(cplR) == 1 :
            """ Dimer to Dimer """
            RG.add_weighted_edges_from([(row, col, pDF[col][row])])
            snodes[row]=1
            snodes[col]=1
          elif len(cplR) == 2 :
            """ 2xMon to 2xMon """
            if cplR[0] != cplC[0] and cplR[1] != cplC[1] : 
              """ this fucks up simulations, skip it! """
              continue
            if cplR[0] != cplC[0] :
              #print "adding", cplR[0], cplC[0]
              edict = RG.get_edge_data(cplR[0], cplC[0])
              if edict :
                if edict[0]['weight'] < pDF[col][row] :
                  edict[0]['weight'] = pDF[col][row]
              else :
                RG.add_weighted_edges_from([(cplR[0], cplC[0], pDF[col][row])])
              snodes[cplR[0]]=1
              snodes[cplC[0]]=1
            if cplR[1] != cplC[1] :
              #print "adding", cplR[1], cplC[1], pDF[col][row], row, col 
              edict = RG.get_edge_data(cplR[1], cplC[1])
              if edict :
                if edict[0]['weight'] < pDF[col][row] :
                  edict[0]['weight'] = pDF[col][row]
              else :
                RG.add_weighted_edges_from([(cplR[1], cplC[1], pDF[col][row])])
              snodes[cplR[1]]=1
              snodes[cplC[1]]=1
          else :
            sys.exit('no support for compexes of size > 2')
        else :
          hcount[0] += 1
          hnode = ('HYPER', str(hcount[0]))
          for tup in cplR+cplC :
            if not RG.has_node(tup):
              #print tup, "recovered!"
              snodes[tup]=1
          #print "Adding:", row, col, pDF[col][row]
          for tup in cplR :
            RG.add_weighted_edges_from([(tup, hnode, pDF[col][row])])
          for tup in cplC :
            RG.add_weighted_edges_from([(hnode, tup, pDF[col][row])])
  else :
    print """ specify valid modus (uni/bi) """
  
  return

def get_complexes((spe, ss)) :
  """ Split a complex of interacting secondary structures into the connected
  components: Input e.g. in form tuple: ('A&B', '.(...).&....')

  :return: a sorted list of connected components [('A',.(...).'),('B','....')]
  """
  species   = zip(spe.split('&'), ss.split('&'))
  complexes = []
  graph     = nx.Graph()

  """ Initialize a disconnected graph including symmetric homodimers! """
  for e, mol in enumerate(species) :
    graph.add_node((e, mol))

  s = 0; stack=[]
  for char in ss:
    molA = (s, species[s])
    if (char == '('):
      stack.append(molA);
    elif (char == ')'):
      molB=stack.pop();
      graph.add_edge(molA,molB)
    elif (char == '&'):
      s += 1

  """ Sorting to access the index 's' in the sets of connected components,
  that is important for homodimers """
  for cocp in sorted(list(nx.connected_components(graph)), 
      key=lambda x : list(x)[0][0]): 
    stru = ''
    comp = ''
    for c,(ch,ss) in sorted(list(cocp), key=lambda x: x[0], reverse=False):
      stru += ss+'&'
      comp += ch+'&'
    stru = stru[:-1]
    comp = comp[:-1]
    complexes.append((comp, stru))

  return sorted(complexes, key=lambda x:x[0])

def get_interkin_args():
  """ A collection of arguments that are used by **interkin.py** """
  parser = argparse.ArgumentParser()

  parser.add_argument("--RNAsubopt", default='RNAsubopt', action = 'store',
      help="Specify path to your *RNAsubopt* executable")
  parser.add_argument("--barriers", default='barriers', action = 'store',
      help="Specify path to your *barriers* executable") 
  parser.add_argument("--treekin", default='treekin', action = 'store',
      help="Specify path to your *treekin* executable")
  parser.add_argument("--tmpdir", default='BarMap_Files', action = 'store',
      help="Specify path for storing temporary output files")

  parser.add_argument("-v","--verbose", action="store_true",
      help="Track process by writing verbose output during calculations") 
  parser.add_argument("-f","--force", action="store_true",
      help="Force to overwrite existing BarMap files") 
  parser.add_argument("-n","--name", default='',
      help="Name your output files, this option overwrites the fasta-header")
  #parser.add_argument("--circ", action="store_true",
  #    help="Circular RNA")  
  #parser.add_argument("--noLP", action="store_true",
  #    help="Compute BarMap pipeline with the ViennaRNA --noLP option")
  parser.add_argument("-T","--temperature", type=float, default=37.0,
      help="Set the temperature for ViennaRNA computations")

  #TODO: need to think whether this is the best way to do it... hidden, for now!
  parser.add_argument("--spatch", 
      default='pylands_spatch.py', action = 'store',
      help=argparse.SUPPRESS)
      #help="Specify a script to postprocess suboptimals")
  #parser.add_argument("--s_stmp", default='/tmp', action = 'store',
  #    help="Specify path to a temporary sort directory for unix sort")

  parser.add_argument("-e", "--s_ener", type=float, default=0,
      help="Set the energy range for suboptimal structure computation." + 
      " This will overwrite the options --s_maxe and --s_maxn.")
  parser.add_argument("--s_maxe", type=float, default=20,
      help="Set a the maximum subopt range in combination with --s_maxn.")
  parser.add_argument("--s_maxn", type=int, default=7100000,
      help="Specify the number of suboptimal structures. The corresponding" +
      " energy range is computed form the full length molecule.")
  parser.add_argument("--b_minh", type=float, default=0.001,
      help="Set the minimum barrier height (i.e. barriers --minh)")
  parser.add_argument("--b_maxn", type=int, default=100,
      help="Set the maximum number of local minima (i.e. barriers --max)")

  #parser.add_argument("--p0", nargs='+', default=['1=1'])
  parser.add_argument("--k0", type=float, default=2e5)
  #parser.add_argument("--t0", type=float, default=0.0)
  #parser.add_argument("--ti", type=float, default=1.02)
  #parser.add_argument("--t8", type=float, default=0.02)
  #parser.add_argument("--tX", type=float, default=60)
  #parser.add_argument("--cutoff", type=float, default=0.01)

  parser.add_argument("--nohomo", 
      help="Turn off homo-dimers", action="store_true")
  parser.add_argument("--d3server", 
      help="Start an interactive viz", action="store_true")
  parser.add_argument("--symmetry", 
      help="correct homo-dimers for their rotational symmetry", 
      action="store_true")
 
  return parser.parse_args()

def main():
  """ Interkin - folding kinetics of two interacting RNA molecules
    * call the pipeline: "RNAsubopt | barriers" in order to get tree files and
        rate matrices
    * make a reaction graph composed of uni- and bi-molecular reactions found
        in the rate matrices
    * visualize the graph using networkx or D3js
    * calculate the equilibrium distribution and see how much it differs from 
        the complete cofolded ensemble
    * dump the graph in order to translate it into a ODE system using the perl
        SundialsWrapper library

    TODO: write a landscape library in order to get rid of BarMap dependency
  """
  args = get_interkin_args()

  """ Read Input & Update Arguments """
  name, inseq = nal.parse_vienna_stdin()

  if args.name == '' : 
    args.name = name
  else :
    name = args.name

  if args.verbose: 
    print "# Input: {:s} {:s} {:6.2f} kcal/mol".format(name, inseq, args.s_ener)
    
  if args.s_ener == 0 :
    args.s_ener, args.s_maxn = nal.sys_subopt_range(inseq, 
        nos=args.s_maxn, 
        maxe=20, #args.s_maxe, 
        verb=args.verbose)
    if args.verbose:
      print "# Energyrange {:.2f} computes {:d} sequences".format(
          args.s_ener, args.s_maxn)

  """ Dirty way to use RNAsubopt symmetry correction as commandline interface """
  s_options = ['|', args.spatch]
  if args.symmetry : s_options.append('-s') 

  """ Choose the set of species for further analysis """
  # default: [A, A, AA], [A, B, AB, AA, BB], [A, B, C, AA, AB, AC, BB, BC, CC] 
  # nohomo: [A, B, AB] (for debugging)
  specieslist = get_species(inseq, args.nohomo)

  """ Compute folding kinetics """
  bionly = False # Debugging option, compute kinetics only from cofold tree
  RG, nodes = cofold_barriers(name, specieslist, 
      s_ener = args.s_ener,
      s_maxn = args.s_maxn,
      s_opts = s_options,
      barriers='barriers',
      minh=args.b_minh,
      maxn=args.b_maxn,
      k0=args.k0,
      temp=37.0,
      moves='single-base-pair',
      gzip=True,
      rates=True,
      bsize=False,
      saddle=False,
      bionly=bionly,
      force=args.force,
      verb=args.verbose)

  dump_graph(RG, nodes, args)

  '''

  """ Compute equilibrium properties """
  # plot RNAcofold equilibrium properties with a given range of a0s and b0s
  # plot gradient basin equilibrium for a given range of a0s and b0s
  # => [seq&seq], [K]
  kT=0.61632077549999997
  a0s = map(lambda x: 10**-x, range(9,8,-3))
  b0s = map(lambda x: 10**-x, range(6,1,-1))

  eA =0
  eB =0
  eAA=0
  eAB=0
  eBB=0
  for (name, species) in specieslist :
    if '&' not in name: continue
    [seqA, seqB] = species.split('&')
    [nA, nB] = name.split('&')
    plotname = args.name+'_'+nA+'_'+nB
    homodimer = (seqA == seqB)

    """ if it is a homodimer, G_A = G_B and G_AB = G_AA = G_BB """
    if not homodimer :
      [G_A, G_B, G_AB, G_AA, G_BB] = conc.cofold_free_energies(species)
      print 'ref_cfold', species, nA, G_A, nB, G_B, \
          nA+nB, G_AB, nA*2, G_AA, nB*2, G_BB
      eqdata = conc.cofold_equilibrium(species, a0s, b0s)
      conc.cofold_plot(eqdata,plotname+'_cofold_plot')

      """
      [G_A, G_B, G_AB, G_AA, G_BB] = basin_free_energies(species)
      print 'ref_cfold', species, nA, G_A, nB, G_B, \
          nA+nB, G_AB, nA*2, G_AA, nB*2, G_BB
      eqdata = conc.cofold_equilibrium(species, a0s, b0s)
      """
    
    RNA.cvar.cut_point = len(seqA)+1
    [stmp,G_A,G_B,G_AB,waste] = RNA.co_pf_fold(seqA+seqB,None)
    K = math.e**((G_A+G_B-G_AB)/kT)
    print 'ref_bimol', species, nA, G_A, nB, G_B, nA+nB, G_AB, "K", K
    if homodimer :
      eqdata = conc.bimol_equilibrium(K, b0s, [])
      conc.homodim_plot(eqdata,plotname+'_bimol_plot')
    else :
      eqdata = conc.bimol_equilibrium(K, a0s, b0s)
      conc.bimol_plot(eqdata,plotname+'_bimol_plot')

    if not bionly :
      [G_A, G_B, G_AB] = unimolecular_basins(RG, name, species, args, nodes)
      K = math.e**((G_A+G_B-G_AB)/kT)
      print 'uitree', species, nA, G_A, nB, G_B, nA+nB, G_AB, "K", K

    [G_ApB, G_AB] = bimolecular_cofold_basins(
        RG, name, species, args, nodes, bionly)
    K = math.e**((G_ApB-G_AB)/kT)
    print 'bitree', species, nA+'+'+nB, G_ApB, nA+nB, G_AB, "K", K
    if homodimer :
      eqdata = conc.bimol_equilibrium(K, b0s, [])
      conc.homodim_plot(eqdata, plotname+'_coarse_grain_plot')
    else :
      eqdata = conc.bimol_equilibrium(K, a0s, b0s)
      [a0, b0, eA, eB, eAB] = zip(*eqdata)
      conc.bimol_plot(eqdata, plotname+'_coarse_grain_plot')
   
  if False : # Depricated and a bit ugly code, it is for debugging only
    """ Compare lists of gradient basin energies for fake and true dimers """
    interkin_equilibrium(specieslist, args)

  # Calculate equilibrium distribution by hand... (like for prions!)
  # use that as a start for the simulation: 
  # *) calculate the concentrations of species
  # *) calculate the occupancy of an lmin at equilibrium
  # *) return the equilibrium starting population for a simulation ...
  '''

  #if False :
  #  plt.clf()
  #  nx.draw_networkx_labels(RG,pos=nx.random(RG))
  #  plt.savefig(args.name+'_nwx.pdf')
  #if args.d3server :
  #  draw_d3_graph(RG)

  return
  
if __name__ == '__main__':
  main()

