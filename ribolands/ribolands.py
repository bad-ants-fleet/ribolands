#
# ribolands.ribolands
# 
# Home of the RiboLandscape object.
#

import logging
rlog = logging.getLogger(__name__)

import RNA 
import math
from struct import pack
from natsort import natsorted
from crnsimulator import ReactionGraph
from crnsimulator.crn_parser import parse_crn_string

from ribolands.syswraps import sys_treekin
from ribolands.utils import Species, tarjans

class RiboLandscapeError(Exception):
    pass

class RiboLandscape:
    """ Implemented for unimolecular reactions.

    The backbone of any ribolands landscape object, inspired by networkx
    interface, although networkx dependency has been removed.

    Suggested attributes for nodes: structure, identity, energy, occupancy, active
    Suggested attributes for edges: weight
    """
    def __init__(self, sequence, vrna_md = None, prefix = 'n'):
        self.sequence = sequence
        if vrna_md:
            self.md = vrna_md
            self.fc = RNA.fold_compound(sequence, vrna_md)
        else:
            self.md = None
            self.fc = None
        self.prefix = prefix  # for autogenerated node IDs
        self.nodeID = 0       # for autogenerated node IDs
        self._nodes = dict()
        self._edges = dict()
        self._in_edges = dict()
        self._out_edges = dict()

    @property
    def RT(self):
        RT = 0.61632077549999997
        if self.md.temperature != 37.0:
            kelvin = 273.15 + self.md.temperature
            RT = (RT / 310.15) * kelvin
        return RT

    @property
    def nodes(self):
        return self._nodes

    @property
    def edges(self):
        return self._edges

    def has_node(self, n):
        return n in self._nodes

    def has_edge(self, s1, s2):
        return (s1, s2) in self._edges

    def predecessors(self, node):
        for n in self._in_edges[node]:
            yield n

    def successors(self, node):
        for n in self._out_edges[node]:
            yield n

    def addnode(self, key, 
                structure = None, 
                identity = None, 
                energy = None, 
                occupancy = 0, 
                **kwargs):
        """ Add a node with attributes to the graph.

        Key and identity are separate, because there are circumstances where
        one wants to have e.g. the secondary structure as key, but that does
        not go well with visualization or simulation. 
        """
        assert key not in self._nodes
        if energy is not None:
            assert isinstance(energy, int)
            energy = energy
        elif structure is not None:
            energy = int(round(self.fc.eval_structure(structure)*100))
        if identity is None:
            identity = f'{self.prefix}{self.nodeID}'
            self.nodeID += 1
        assert isinstance(identity, str)
        self._nodes[key] = {'structure': structure,
                            'identity': identity,
                            'energy': energy,
                            'occupancy': occupancy}
        self._nodes[key].update(kwargs)
        self._in_edges[key] = set()
        self._out_edges[key] = set()
        return

    def addedge(self, n1, n2, weight = None, **kwargs):
        assert n1 in self._nodes
        assert n2 in self._nodes
        if (n1, n2) not in self._edges:
            self._edges[(n1, n2)] = {'weight': weight}
            self._out_edges[n1].add(n2)
            self._in_edges[n2].add(n1)
        self._edges[(n1, n2)].update(kwargs)
        return

    def sorted_nodes(self, attribute = 'energy', rev = False, nodes = None):
        """ Provide active nodes or new nodes, etc. if needed. """
        if nodes is None:
            nodes = self.nodes
        return sorted(nodes, key = lambda x: self.nodes[x][attribute], reverse = rev)

    def sccs(self, minh = None):
        def neighbors(node):
            for nbr in self.successors(node):
                if minh is not None:
                    se = self.edges[(node, nbr)]['saddle_energy']
                    e1 = self.nodes[node]['energy']
                    if se - e1 > minh:
                        continue
                yield nbr
        species = {n: Species(n) for n in self.nodes}
        products = {v: list(species[vv] for vv in neighbors(k)) for k, v in species.items()}
        for scc in tarjans(list(species.values()), products):
            yield [n.name for n in scc]

    def get_rate(self, s1, s2):
        """ Returns the direct transition rate of two secondary structures. """
        return self._edges[(s1, s2)]['weight'] if self.has_edge(s1, s2) else 0

    def to_crn(self, filename = None):
        def crn_gen():
            for (x, y), data in self.edges.items():
                assert not isinstance(self.nodes[x]['identity'], int)
                assert not isinstance(self.nodes[y]['identity'], int)
                reactant = self.nodes[x]['identity']
                product = self.nodes[y]['identity']
                if data['weight'] == 0:
                    continue
                yield "{:s} -> {:s} [k = {:g}]\n".format(reactant, product, data['weight'])
        if filename:
            with open(filename, 'w') as crn:
                for x in crn_gen():
                    crn.write(x)
            return
        return "".join(crn_gen())

    def to_crnsimulator(self, filename, sorted_vars = None):
        """
        """
        if filename[-3:] != '.py':
            filename += '.py' 
        def irrev(rev):
            # Split CRN into irreversible reactions
            new = []
            for [r, p, k] in rev:
                if None in k:
                    rlog.warning('# Set missing rates to 1.')
                    k[:] = [x if x is not None else 1 for x in k]
                if len(k) == 2:
                    new.append([r, p, k[0]])
                    new.append([p, r, k[1]])
                else:
                    new.append([r, p, k[0]])
            return new

        crnstring = self.to_crn()
        crn, species = parse_crn_string(crnstring)
        RG = ReactionGraph(irrev(crn))
        V = sorted_vars if sorted_vars else natsorted(species)

        # ********************* #
        # PRINT ODE TO TEMPLATE #
        # ..................... #
        filename, _ = RG.write_ODE_lib(sorted_vars = V, concvect = None,
                                             jacobian = False,
                                             filename = filename)
        rlog.info(f'# Wrote ODE system: {filename}')
        return filename

    def simulate_crn(self, filename = None, sorted_vars = None):
        # Take odeint parameters
        # Take plotting parameters
        # call to_crnsimulator and do the I/O handling
        raise NotImplementedError

    def get_simulation_files_tkn(self, basename, snodes = None):
        """ Print a rate matrix and the initial occupancy vector.

        This function prints files and parameters to simulate dynamics using the
        commandline tool treekin. A *.bar file contains a sorted list of present
        structures, their energy and their neighborhood and the corresponding
        energy barriers. A *.rts or *.rts.bin file contains the matrix of
        transition rates either in text or binary format. Additionaly, it returns
        a vector "p0", which contains the present occupancy of structures. The
        order or elements in p0 contains

        Note:
          This routine contains a hack for *.bar files to contain the energy
          barriers to all other nodes in the graph. This is only enabled when
          a get_saddle() attribute is defined for the object (see PrimePathLandscape).

        Args:
            basename (str): Basename of output files.

        Returns:
            [str, str, str, str]: Binary rates file, Text rates file, barriers-like output, p0 string
        """
        seq = self.sequence
        if snodes is None:
            snodes = self.sorted_nodes(attribute = 'energy')
        num = len(snodes) + 1

        bofile = basename + '_lands.bar'
        brfile = basename + '_rates.txt'
        bbfile = basename + '_rates.bin'
        p0 = []

        with open(bofile, 'w') as bar, open(brfile, 'w') as rts, open(bbfile, 'wb') as brts:
            barhead = ' '.join(map("{:7d}".format, range(1, num))) if hasattr(self, 'get_saddle') else ''
            bar.write("  ID {}  Energy  {}\n".format(seq, barhead))
            brts.write(pack("i", len(snodes)))
            for ni, node in enumerate(snodes, 1):
                ns = self.nodes[node]['structure']
                ne = self.nodes[node]['energy']/100
                no = self.nodes[node]['occupancy']
                
                # Calculate barrier heights to all other basins.
                barstr = ''
                if hasattr(self, 'get_saddle'): 
                    for other in snodes:
                        os = self.nodes[other]['structure']
                        oe = self.nodes[other]['energy']/100
                        se = self.get_saddle(node, other)
                        if se is not None:
                            barstr += ' {:7.2f}'.format(se - ne)
                        else:
                            barstr += ' {:7.2f}'.format(float('nan'))

                # Print structures and neighbors to bfile:
                bar.write("{:4d} {} {:7.2f} {}\n".format(ni, ns, ne, barstr))

                # Add ni occupancy to p0
                if no > 0:
                    p0.append("{}={}".format(ni, no))

                # Print rate matrix to rfile and brfile
                trates = []
                rates = []
                for other in snodes:
                    if self.has_edge(node, other):
                        rates.append(self.edges[(node, other)]['weight'])
                    else:
                        rates.append(0)

                    if self.has_edge(other, node):
                        trates.append(self.edges[(other, node)]['weight'])
                    else:
                        trates.append(0)
                line = " ".join(map("{:10.4g}".format, rates))
                rts.write("{}\n".format(line))
                for r in trates:
                    brts.write(pack("d", r))
        return bbfile, brfile, bofile, p0

    def simulate_tkn(self, basename, snodes = None, form = 'pdf', **kwargs):
        bbf, brf, bof, p0 = self.get_simulation_files_tkn(basename, snodes = snodes)
        tfile, _ = sys_treekin(basename, bbf, p0 = p0, **kwargs)
        return tfile

    def __repr__(self):
        return f"{self.__class__.__name__}({self.sequence}, {len(self.nodes)=}, {len(self.edges)=})"

