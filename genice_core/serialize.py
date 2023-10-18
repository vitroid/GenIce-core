"""
Deals with serialization/deserialization of a directed graph.
"""
import networkx as nx


def serialize(dg: nx.DiGraph) -> str:
    """Encode a digraph

    Args:
        dg (nx.DiGraph): the directed graph to be encoded.

    Returns:
        str: a string representing the digraph.
    """
    def encode1(v):
        return "+".join(f"{x}" for x in sorted(dg.successors(v)))

    return ".".join(encode1(v) for v in sorted(dg.nodes()))


def deserialize(s:str) -> nx.Graph:
    """Decode a string into a digraph

    Args:
        s (str): encoded string of a graph

    Returns:
        nx.Graph: the recovered digraph.
    """
    dg = nx.DiGraph()

    for i, neis in enumerate(s.split(".")):
        if neis == "":
            continue
        for j in neis.split("+"):
            dg.add_edge(i, int(j))
    return dg
