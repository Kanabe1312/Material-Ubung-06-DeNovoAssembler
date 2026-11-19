#! /usr/bin/python3
###########################################################################
# Github-Name: Kanbe1312
# Name: Marcu Alexandru Stefan
# Email: alexmarcu670@gmail.com  (github),
# Punkte: xx/17
#############################################################################

class Read:
    def __init__(self, lines):
        self.name = lines[0].strip()[1:]
        self.bases = "".join([x.strip() for x in lines[1:]]).upper()

    def get_kmers(self, kmersize):
        res = {}
        for pos in range(0, len(self.bases) - kmersize + 1):
            kmer = self.bases[pos:(pos + kmersize)]
            if kmer not in res:
                res[kmer] = 0
            res[kmer] += 1
        return res

    def __str__(self):
        return self.name + ": " + self.bases[:20] + "..."

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, Read):
            return False
        if other.name == self.name and other.bases == self.bases:
            return True
        return False


class DBGnode:
    def __init__(self, seq):
        self.seq = seq
        self.eto = dict()
        self.efrom = dict()

    def add_edge_to(self, eto):
        if eto not in self.eto:
            self.eto[eto] = 0
        self.eto[eto] += 1

    def add_edge_from(self, efrom):
        if efrom not in self.efrom:
            self.efrom[efrom] = 0
        self.efrom[efrom] += 1

    def get_potential_from(self):
        res = []
        for x in "AGTC":
            res += [x + self.seq[:-1]]
        return res

    def get_potential_to(self):
        res = []
        for x in "AGTC":
            res += [self.seq[1:] + x]
        return res

    def get_edge_to_weight(self, other):
        if not other in self.eto:
            return 0
        return self.eto[other]

    def get_edge_from_weight(self, other):
        if not other in self.efrom:
            return 0
        return self.efrom[other]


class DBGraph:
    # 2 Punkte
    def __init__(self):
        self.nodes = {}
        self.kmer_size = None
    # 8 Punkte
       # 8 Punkte
    def add_kmers(self, kmers):
        if not kmers:
            return

        for kmer in kmers:
            if self.kmer_size is None:
                self.kmer_size = len(kmer)
            elif len(kmer) != self.kmer_size:
                raise ValueError("Different k-mer sizes")

            if kmer not in self.nodes:
                self.nodes[kmer] = DBGnode(kmer)

        for kmer in kmers:
            node = self.nodes[kmer]

            for potential in node.get_potential_to():
                if potential in self.nodes:
                    target = self.nodes[potential]

                    if target not in node.eto:
                        node.add_edge_to(target)
                        target.add_edge_from(node)

            for potential in node.get_potential_from():
                if potential in self.nodes:
                    source = self.nodes[potential]

                    if source not in node.efrom:
                        node.add_edge_from(source)
                        source.add_edge_to(node)
               
    # 2 Punkte
    def count_edges(self):
        count = 0

        for node in self.nodes.values():
            count += len(node.eto)
        return count    
    # 1 Punkt
    def count_nodes(self):
        return len(self.nodes)
    # 1 Punkt
    def __str__(self):
       return (f"{self.kmer_size} kmer größe\n"
               f"{self.count_nodes()} k-mere/Knoten\n"
               f"{self.count_edges()} Kanten")


def read_fasta(readfile):
    with open( readfile, 'r') as readfile:
        readList = []
        buffer = []
        for line in readfile:
            if line.startswith(">"):
                if buffer:
                    readList.append(Read(buffer))
                buffer = [line]
            else:
                buffer.append(line.strip())
        if buffer:
            readList.append(Read(buffer))
        return readList

# 3 Punkte
def build_graph(filename, kmersize):
    graph = DBGraph()

    reads = read_fasta(filename)

    for read in reads:
        kmers = read.get_kmers(kmersize)
        graph.add_kmers(kmers)

    return graph


if __name__ == "__main__":
    pass
