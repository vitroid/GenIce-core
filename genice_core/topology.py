"""
Arrange edges appropriately.
"""

from logging import getLogger
from typing import Union
import networkx as nx
import numpy as np

from genice_core.dipole import minimize_net_dipole

# Generate documents for these functions only.
__all__ = []


def extend_path(g: nx.Graph, path: list) -> list:
    """extend the path

    Args:
        g (nx.Graph): a linear or a simple cyclic graph.
        path (list): A given path tho be extended

    Returns:
        list: the extended path or cycle
    """
    while True:
        last, head = path[-2:]
        for next in g.neighbors(head):
            if next != last:
                break
        else:
            # no next node
            return path
        path.append(next)
        if next == path[0]:
            # is cyclic
            return path


def find_path(g: nx.Graph) -> list:
    """Find the path in g

    Args:
        g (nx.Graph): a linear or a simple cyclic graph.

    Returns:
        list: the path or cycle
    """
    nodes = list(g.nodes())
    head = nodes[0]
    nei = [v for v in g.neighbors(head)]
    if len(nei) == 0:
        # lone node
        return []
    elif len(nei) == 1:
        # an end node, fortunately.
        return extend_path(g, [head, nei[0]])
    c0 = extend_path(g, [head, nei[0]])

    if c0[-1] != head:
        # linear chain graph
        # join the another side
        c1 = extend_path(g, [head, nei[1]])
        return list(reversed(c0)) + c1[1:]

    # cyclic graph
    return c0


def noodlize(g: nx.Graph) -> nx.Graph:
    """Divide each vertex of the graph and make a set of paths.

    A new algorithm suggested by Prof. Sakuma, Yamagata University.

    Args:
        g (nx.Graph): An ice-like undirected graph. All vertices must not be >4-degree.

    Returns:
        nx.Graph: A graph mode of chains and cycles.
    """

    nnode = len(g)
    # 1. Split the nodes.

    # divided graph
    divg = nx.Graph(g)

    for v in g:
        nei = list(divg.neighbors(v))
        assert len(nei) <= 4, "degree must be <=4"
        # fill by Nones if number of neighbors is less than 4
        nei = (nei + [None, None, None, None])[:4]
        # two neighbor nodes that are passed away to the new node
        migrate = set(np.random.choice(nei, 2, replace=False)) - set([None])
        # new node label
        newv = v + nnode
        # assemble edges
        for x in migrate:
            divg.remove_edge(x, v)
            divg.add_edge(newv, x)

    # divg is made of chains and cycles.
    return divg


def decompose_complex_path(path: list):
    """A generator that divides a complex path with self-crossings to set of simple cycles and paths.

    Args:
        path (list): A complex path

    Yields:
        list: a short and simple path/cycle
    """
    logger = getLogger()
    if len(path) == 0:
        return
    logger.debug(f"decomposing {path}...")
    order = dict()
    order[path[0]] = 0
    store = [path[0]]
    headp = 1
    while headp < len(path):
        node = path[headp]
        if node in order:
            # it is a cycle!
            size = len(order) - order[node]
            cycle = store[-size:] + [node]
            yield cycle
            # remove them from the order[]
            for v in cycle[1:]:
                del order[v]
            # truncate the store
            store = store[:-size]
        order[node] = len(order)
        store.append(node)
        headp += 1
        logger.debug([order, store])
    if len(store) > 1:
        yield store
    logger.debug(f"Done decomposition.")


def simple_paths(
    nnode: int,
    divg: nx.Graph,
):
    """
    Set the orientations to the components.

    nnode (int): number of nodes in the original graph.
    divg (nx.Graph): the divided graph.
    """
    for c in nx.connected_components(divg):
        # a component of c is either a chain or a cycle.
        subg = divg.subgraph(c)
        nn = len(subg)
        ne = len([e for e in subg.edges()])
        assert nn == ne or nn == ne + 1
        # Find a simple path in the doubled graph
        # It must be a simple path or a simple cycle.
        path = find_path(subg)
        # Flatten then path. It may make the path self-crossing.
        path = [v % nnode for v in path]
        # Divide a long path into simple paths and cycles.
        yield from decompose_complex_path(path)
