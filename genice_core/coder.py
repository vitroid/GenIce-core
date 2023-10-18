import networkx as nx

from genice_core.topology import find_path


def encode(dg):
    """Encode the digraph"""

    def encode1(v):
        return "+".join(f"{x}" for x in sorted(dg.successors(v)))

    return ".".join(encode1(v) for v in sorted(dg.nodes()))


def decode(s):
    dg = nx.DiGraph()

    for i, neis in enumerate(s.split(".")):
        if neis == "":
            continue
        for j in neis.split("+"):
            dg.add_edge(i, int(j))
    return dg
