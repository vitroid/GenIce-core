# ノードを分割し、成分に分け、パスをさがし、という処理が無駄に遅い気がするので、自作してみる。

from collections import defaultdict
from logging import getLogger
from typing import Dict, Iterable, MutableSequence, Union
import networkx as nx


class Graph(Dict):
    def __init__(self, g=None):
        # logger = getLogger()
        self._edges = defaultdict(set)
        # logger.info(type(g))
        if type(g) is list:
            self.add_edges(g)
        elif g is None:
            pass
        elif type(g) is Graph:
            self.add_edges(g.edges())
        elif type(g) is DiGraph:
            self.add_edges(g.edges())
        elif type(g) is nx.Graph:
            self.add_edges(g.edges())
        elif type(g) is nx.DiGraph:
            self.add_edges(g.edges())
        else:
            assert False

    def add_edges(self, L):
        for i, j in L:
            self._edges[i].add(j)
            self._edges[j].add(i)

    def add_edge(self, i, j):
        self._edges[i].add(j)
        self._edges[j].add(i)

    def __iter__(self):
        yield from self._edges

    def __getitem__(self, n):
        return self._edges[n]

    # def __setitem__(self):
    #     pass

    def __len__(self):
        return len(self._edges)

    def has_node(self, x):
        return x in self._edges

    def degree(self, x):
        return len(self._edges[x])

    def nodes(self):
        return self._edges.keys()

    def edges(self):
        for key, values in self._edges.items():
            for v in values:
                if key < v:
                    yield key, v

    def neighbors(self, x):
        return self._edges.get(x, set())

    def remove_edge(self, x, y):
        self._edges[x].remove(y)
        self._edges[y].remove(x)

    def subgraph(self, x):
        g = type(self)()  # empty clone
        for i in x:
            for j in self._edges[i]:
                if j in x:
                    g.add_edge(i, j)
        return g


def connected_components(g):
    def mygroup(id):
        while id >= 0:
            id = group[id]
        return -id - 1

    group = {x: -x - 1 for x in g}

    members = {x: set([x]) for x in g}

    for i, j in g.edges():
        hi = mygroup(i)
        hj = mygroup(j)
        if hi != hj:
            members[hi] |= members[hj]
            group[hj] = hi

    for i, value in members.items():
        if group[i] < 0:
            yield value


def add_path(g, L):
    for i, j in zip(L, L[1:]):
        g.add_edge(i, j)


class DiGraph(Graph):
    def add_edges(self, L):
        for i, j in L:
            self._edges[i].add(j)

    def remove_edge(self, x, y):
        self._edges[x].remove(y)

    def add_edge(self, i, j):
        self._edges[i].add(j)

    def edges(self):
        for key, values in self._edges.items():
            for v in values:
                yield key, v

    def nodes(self):
        assert False


if __name__ == "__main__":
    g = Graph([(1, 2), (2, 3), (4, 5)])
    print(g[2])
    for c in connected_components(g):
        print(c)
