
import RNA 
import networkx as nx
from struct import pack

class RiboLandscape(nx.DiGraph):
    """ Implemented for unimolecular reactions.

    Suggested attributes for nodes: energy, structure, identity.
    Suggested attributes for edges: weight, saddleE.
    """

    def __init__(self, sequence, vrna_md):
        super(RiboLandscape, self).__init__()
        self.sequence = sequence
        self.md = vrna_md
        self.fc = RNA.fold_compound(sequence, vrna_md)

    def sorted_nodes(self, attribute = 'energy', rev = False):
        return sorted(self.nodes(), key = lambda x: self.nodes[x][attribute], reverse = rev)

    def get_saddle(self, s1, s2):
        """Returns the saddle energy of a transition edge."""
        if self.has_edge(s1, s2):
            return self[s1][s2]['saddleE']
        else:
            return None

    def get_barrier(self, s1, s2):
        """Returns the barrier energy of a transition edge."""
        if self.has_edge(s1, s2):
            return self[s1][s2]['saddleE'] - self.node[s1]['energy']
        else:
            return None

    def get_rate(self, s1, s2):
        """Returns the direct transition rate of two secondary structures."""
        if self.has_edge(s1, s2):
            return self[s1][s2]['weight']
        else:
            return 0

    def get_simulation_files_tkn(self, basename):
        """ Print a rate matrix and the initial occupancy vector.

        This function prints files and parameters to simulate dynamics using the
        commandline tool treekin. A *.bar file contains a sorted list of present
        structures, their energy and their neighborhood and the corresponding
        energy barriers. A *.rts or *.rts.bin file contains the matrix of
        transition rates either in text or binary format. Additionaly, it returns
        a vector "p0", which contains the present occupancy of structures. The
        order or elements in p0 contains

        Note:
          A *.bar file contains the energy barriers to all other nodes in the
          graph. So it is not the same as a "classic" barfile produce by
          barriers.

        Args:
            basename (str): Basename of output files.

        Returns:
            [str, str, str]: Binary rates file, Text rates file, barriers-like output.
        """
        seq = self.sequence

        snodes = self.sorted_nodes(attribute = 'energy')
        num = len(snodes) + 1

        bofile = basename + '_lands.bar'
        brfile = basename + '_rates.txt'
        bbfile = basename + '_rates.bin'

        with open(bofile, 'w') as bar, open(brfile, 'w') as rts, open(bbfile, 'wb') as brts:
            bar.write("  ID {}  Energy  {}\n".format(seq, 
                ' '.join(map("{:7d}".format, range(1, num)))))
            brts.write(pack("i", len(snodes)))
            for ni, node in enumerate(snodes):
                ns = self.nodes[node]['structure']
                ne = self.nodes[node]['energy']
                
                # Calculate barrier heights to all other basins.
                barstr = ''
                for other in snodes:
                    os = self.nodes[other]['structure']
                    oe = self.nodes[other]['energy']
                    sE = self.get_saddle(node, other)
                    if sE is not None:
                        barstr += ' {:7.2f}'.format(sE - ne)
                    else:
                        barstr += ' {:7.2f}'.format(float('nan'))

                # Print structures and neighbors to bfile:
                bar.write("{:4d} {} {:7.2f} {}\n".format(ni, ns, ne, barstr))

                # Print rate matrix to rfile and brfile
                trates = []
                rates = []
                for other in snodes:
                    if self.has_edge(node, other):
                        rates.append(self[node][other]['weight'])
                        trates.append(self[other][node]['weight'])
                    else:
                        rates.append(0)
                        trates.append(0)
                line = " ".join(map("{:10.4g}".format, rates))
                rts.write("{}\n".format(line))
                for r in trates:
                    brts.write(pack("d", r))

        return bbfile, brfile, bofile

    def to_crn(self, basename):
        # An easier way to simulate it ...
        pass

    def get_simulator(self, basename):
        # simple version for unimolecular reactions only? 
        # or shall we support the (k1, k2) model?
        pass

