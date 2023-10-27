"""
Arrange edges appropriately.
"""

from logging import getLogger, DEBUG
import networkx as nx
import numpy as np
from typing import Union
import random

# Generate documents for these functions only.
# __all__ = []


def _extend_path(g: nx.Graph, path: list) -> list:
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


def _find_path(g: nx.Graph) -> list:
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
        return _extend_path(g, [head, nei[0]])
    c0 = _extend_path(g, [head, nei[0]])

    if c0[-1] != head:
        # linear chain graph
        # join the another side
        c1 = _extend_path(g, [head, nei[1]])
        return list(reversed(c0)) + c1[1:]

    # cyclic graph
    return c0


def _divide(divg: nx.Graph, nei: list, vertex: int, offset: int):
    # fill by Nones if number of neighbors is less than 4
    nei = (nei + [None, None, None, None])[:4]
    # two neighbor nodes that are passed away to the new node
    migrants = set(np.random.choice(nei, 2, replace=False)) - set([None])
    # new node label
    newVertex = vertex + offset
    # assemble edges
    for migrant in migrants:
        divg.remove_edge(migrant, vertex)
        divg.add_edge(newVertex, migrant)


def _divide_node(
    divg: nx.Graph, vertex: int, offset: int, fixed: nx.DiGraph, numFixedEdges: int
):
    nei = list(divg.neighbors(vertex))
    fullDegree = len(nei) + numFixedEdges
    # minimal assertion
    assert fullDegree <= 4, f"degree {fullDegree} > 4"

    # 一部の辺が固定されている場合、分割は自由に行えない。この部分の処理が非常にこみいったことになる。
    if numFixedEdges == 0:
        _divide(divg, nei, vertex, offset)
    elif numFixedEdges == 1:
        # 自由結合が3本: 1本が固定されているので、一本が分割後に固定されることになる。
        # 本来なら前処理でそのような状況をなくしておかれるべき。
        assert len(nei) != 3

        # 前処理のあとで自由結合が2本の場合、第4のダングリング結合が固定結合をキャンセルしているはずなので、分割不要。

        # 自由結合が1の場合はそもそも分割不要。
    elif numFixedEdges == 2:
        # 2本の固定結合が1i1oの場合のみ、残る2本は自由に向きを決められる。
        # 2i/2oの場合には、前処理で残る2本も固定されていなければならないので、ここにはこない。
        pass
        # 残りは割らなくていい。
    elif numFixedEdges == 3:
        # 残る1本は自動的に向きが定まるので前処理で排除されていなければならない。
        # 可能なのは3結合で3本固定の場合のみ。
        assert len(nei) == 0


def noodlize(g: nx.Graph, fixed: Union[nx.DiGraph, None] = nx.DiGraph()) -> nx.Graph:
    """Divide each vertex of the graph and make a set of paths.

    A new algorithm suggested by Prof. Sakuma, Yamagata University.

    Args:
        g (nx.Graph): An ice-like undirected graph. All vertices must not be >4-degree.
        fixed (Union[nx.DiGraph,None], optional): Specifies the edges whose direction is fixed.. Defaults to None.

    Returns:
        nx.Graph: A graph mode of chains and cycles.
    """

    logger = getLogger()

    fixg = nx.Graph(fixed)  # undirected copy

    offset = len(g)

    # divided graph
    divg = nx.Graph(g)
    for edge in fixed.edges():
        divg.remove_edge(*edge)

    for v in g:
        if fixg.has_node(v):
            nfixed = fixg.degree[v]
        else:
            nfixed = 0
        _divide_node(divg, v, offset, fixed, nfixed)

    # divg is made of chains and cycles.
    # divg does not contain the edges in fixed.
    return divg


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
    divg: nx.Graph,
):
    """Set the orientations to the components.

    Args:
        nnode (int): number of nodes in the original graph.
        divg (nx.Graph): the divided graph.

    Yields:
        list: a short and simple path/cycle
    """

    for c in nx.connected_components(divg):
        # a component of c is either a chain or a cycle.
        subg = divg.subgraph(c)
        nn = len(subg)
        ne = len([e for e in subg.edges()])
        assert nn == ne or nn == ne + 1
        # Find a simple path in the doubled graph
        # It must be a simple path or a simple cycle.
        path = _find_path(subg)
        # Flatten then path. It may make the path self-crossing.
        path = [v % nnode for v in path]
        # Divide a long path into simple paths and cycles.
        yield from _decompose_complex_path(path)


# def _decorate(g: nx.Graph):
#     """Add virtual nodes to make all nodes 4-connected.

#     Args:
#         g (nx.Graph): _description_
#     """
#     nodes = list(g)
#     for i in nodes:
#         room = 4 - len(g[i])
#         for j in range(room):
#             g.add_edge(i, -j - 1)


def _remove_dummy_nodes(g: Union[nx.Graph, nx.DiGraph]):
    for i in range(-1, -5, -1):
        if g.has_node(i):
            g.remove_node(i)


def balance(fixed: nx.DiGraph, g: nx.Graph, hook=None):
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
        # draw(pos, g, fixed)
        # 水分子は仮想的に4本の水素結合を持っている。
        # Balance principleにより、
        # 入結合が固定されたら、出結合も同じ数だけ固定されねばならない。
        neighborNodes = list(g[node]) + [-1, -2, -3, -4]
        neighborNodes = neighborNodes[:4]
        while fixed.in_degree(node) > fixed.out_degree(node):
            next = random.choice(neighborNodes)
            if not (fixed.has_edge(node, next) or fixed.has_edge(next, node)):
                fixed.add_edge(node, next)
                perimeter.append(next)
        while fixed.in_degree(node) < fixed.out_degree(node):
            next = random.choice(neighborNodes)
            if not (fixed.has_edge(node, next) or fixed.has_edge(next, node)):
                fixed.add_edge(next, node)
                perimeter.append(next)

    # これで残った無向グラフにgenice-coreして、ちゃんとice ruleが満足されることを示す。

    # ダミーノードを除去する
    # undecorate(g)
    _remove_dummy_nodes(fixed)
