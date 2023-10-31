"""
Arrange edges appropriately.
"""

from logging import getLogger, DEBUG
import networkx as nx
import numpy as np
from typing import Union


def _trace_path(g: nx.Graph, path: list) -> list:
    """trace the path

    Args:
        g (nx.Graph): a linear or a simple cyclic graph.
        path (list): A given path tho be extended

    Returns:
        list: the extended path or cycle
    """
    while True:
        # look at the head of the path
        last, head = path[-2:]
        for next in g[head]:
            if next != last:
                # go ahead
                break
        else:
            # no next node
            return path
        path.append(next)
        if next == path[0]:
            # is cyclic
            return path


def _find_path(g: nx.Graph) -> list:
    """Find the path in g

    Args:
        g (nx.Graph): a linear or a simple cyclic graph.

    Returns:
        list: the path or cycle
    """
    nodes = list(g.nodes())
    # choose one node
    head = nodes[0]
    # look neighbors
    neighbors = list(g[head])
    if len(neighbors) == 0:
        # isolated node
        return []
    elif len(neighbors) == 1:
        # head is an end node, fortunately.
        return _trace_path(g, [head, neighbors[0]])
    # look forward
    c0 = _trace_path(g, [head, neighbors[0]])

    if c0[-1] == head:
        # cyclic graph
        return c0

    # look backward
    c1 = _trace_path(g, [head, neighbors[1]])
    return c0[::-1] + c1[1:]


def _divide(g: nx.Graph, vertex: int, offset: int):
    # fill by Nones if number of neighbors is less than 4
    nei = (list(g[vertex]) + [None, None, None, None])[:4]

    # two neighbor nodes that are passed away to the new node
    migrants = set(np.random.choice(nei, 2, replace=False)) - set([None])

    # new node label
    newVertex = vertex + offset

    # assemble edges
    for migrant in migrants:
        g.remove_edge(migrant, vertex)
        g.add_edge(newVertex, migrant)


def noodlize(g: nx.Graph, fixed: nx.DiGraph = nx.DiGraph()) -> nx.Graph:
    """Divide each vertex of the graph and make a set of paths.

    A new algorithm suggested by Prof. Sakuma, Yamagata University.

    Args:
        g (nx.Graph): An ice-like undirected graph. All vertices must not be >4-degree.
        fixed (Union[nx.DiGraph,None], optional): Specifies the edges whose direction is fixed.. Defaults to None.

    Returns:
        nx.Graph: A graph mode of chains and cycles.
    """

    logger = getLogger()

    g_fix = nx.Graph(fixed)  # undirected copy

    offset = len(g)

    # divided graph
    g_noodles = nx.Graph(g)
    for edge in fixed.edges():
        g_noodles.remove_edge(*edge)

    for v in g:
        if g_fix.has_node(v):
            nfixed = g_fix.degree[v]
        else:
            nfixed = 0
        if nfixed == 0:
            _divide(g_noodles, v, offset)

    # divg is made of chains and cycles.
    # divg does not contain the edges in fixed.
    return g_noodles


def _decompose_complex_path(path: list):
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


def split_into_simple_paths(
    nnode: int,
    g_noodles: nx.Graph,
):
    """Set the orientations to the components.

    Args:
        nnode (int): number of nodes in the original graph.
        divg (nx.Graph): the divided graph.

    Yields:
        list: a short and simple path/cycle
    """

    for verticeSet in nx.connected_components(g_noodles):
        # a component of c is either a chain or a cycle.
        g_noodle = g_noodles.subgraph(verticeSet)
        # nn = len(g_noodle)
        # ne = len([e for e in g_noodle.edges()])
        # assert nn == ne or nn == ne + 1

        # Find a simple path in the doubled graph
        # It must be a simple path or a simple cycle.
        path = _find_path(g_noodle)

        # Flatten then path. It may make the path self-crossing.
        flatten = [v % nnode for v in path]

        # Divide a long path into simple paths and cycles.
        yield from _decompose_complex_path(flatten)


def _remove_dummy_nodes(g: Union[nx.Graph, nx.DiGraph]):
    for i in range(-1, -5, -1):
        if g.has_node(i):
            g.remove_node(i)


def balance(fixed: nx.DiGraph, g: nx.Graph, hook=None):
    """Extend the prefixed digraph to make the remaining graph balanced.

    Args:
        fixed (nx.DiGraph): fixed edges
        g (nx.Graph): skeletal graph
    """
    # prepare the perimeter
    perimeter = [
        node
        for node in fixed
        if fixed.in_degree[node] + fixed.out_degree[node] < g.degree[node]
    ]

    while len(perimeter) > 0:
        node = perimeter.pop(0)
        if node < 0:
            continue
        if hook is not None:
            hook(fixed)

        # fill if degree is less than 4
        neighborNodes = list(g[node]) + [-1, -2, -3, -4]
        neighborNodes = neighborNodes[:4]

        while fixed.in_degree(node) > fixed.out_degree(node):
            next = np.random.choice(neighborNodes)
            if not (fixed.has_edge(node, next) or fixed.has_edge(next, node)):
                fixed.add_edge(node, next)
                perimeter.append(next)

        while fixed.in_degree(node) < fixed.out_degree(node):
            next = np.random.choice(neighborNodes)
            if not (fixed.has_edge(node, next) or fixed.has_edge(next, node)):
                fixed.add_edge(next, node)
                perimeter.append(next)

    _remove_dummy_nodes(fixed)
