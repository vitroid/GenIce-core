"""
固定辺の処理。

とりあえず試験。
"""
import pairlist as pl
import numpy as np
import networkx as nx
import yaplotlib as yap
import random
from logging import getLogger, basicConfig, INFO, DEBUG
import genice_core


def draw(pos, g, dg):
    s = yap.Size(0.2)
    s += yap.ArrowType(2)
    for i in g:
        if i >= 0:
            for j in g[i]:
                if j >= 0:
                    s += yap.Color(0)  # black
                    s += yap.Line(pos[j], pos[i])
    for i in dg:
        if i >= 0:
            for j in dg[i]:
                if j >= 0:
                    s += yap.Color(3)  # green
                    s += yap.Arrow(pos[j], pos[i])
    print(s)  # with new page


def diamond(N: int) -> np.ndarray:
    """Diamond lattice.

    Args:
        N (int): Number of unit cells per an edge of the simulation cell.

    Returns:
        np.ndarray: _description_
    """
    # make an FCC
    xyz = (
        np.array(
            [
                (x, y, z)
                for x in range(N)
                for y in range(N)
                for z in range(N)
                if (x + y + z) % 2 == 0
            ]
        )
        * 2
    )
    xyz = np.vstack([xyz, xyz + 1])
    return xyz


logger = getLogger()
basicConfig(
    level=INFO,
    style="{",
    format="{levelname}:{filename}:{funcName}:{lineno}:{message}",
)


# random.seed(0)

# diamond graph
N = 8
pos = diamond(N)

# 丸く切りとる
center = np.array([N, N, N]) + 0.3
pos = pos[np.linalg.norm(pos - center, axis=1) < N + 0.1]

# 無向グラフ。
g = nx.Graph(
    [
        (i, j)
        for i, j in pl.pairs_iter(pos, 3**0.5 + 0.1, fractional=False, distance=False)
    ]
)

# 結合が4本に満たない場合は負数を入れておく。
# 仮想的に、すべての水分子が4配位であるものとする。
# (あとの処理が簡単にできるので)
# でも表面積が大きいと時間がばかにならない。なしにできないかな。
# できそうだね。

genice_core.topology.decorate(g)

# 固定辺のグラフ。こちらも仮想隣接分子を置く。
fixed = nx.DiGraph()

NH4, F = random.sample(range(len(pos)), 2)

neis = list(g[NH4])
for nei in neis:
    fixed.add_edge(NH4, nei)

neis = list(g[F])
for nei in neis:
    fixed.add_edge(nei, F)


genice_core.topology.balance(fixed, g, hook=lambda fixed: draw(pos, g, fixed))


totalPol = genice_core.dipole.vector_sum(fixed, pos)
logger.info(f"{totalPol}={np.linalg.norm(totalPol):.3f} net dipole before optimization")

# 固定されなかった部分を作る。固定部分の分極をキャンセルするように。
dg = genice_core.ice_graph(
    g,
    vertexPositions=pos,
    dipoleOptimizationCycles=100,
    fixed=fixed,
)
draw(pos, nx.Graph(), dg)

for i, j in fixed.edges:
    assert i >= 0 and j >= 0, (i, j)
    dg.add_edge(i, j)


totalPol = genice_core.dipole.vector_sum(dg, pos)
logger.info(f"{totalPol}={np.linalg.norm(totalPol):.3f} net dipole after optimization")


draw(pos, nx.Graph(), dg)


for node in dg:
    if node not in (NH4, F):
        assert dg.in_degree[node] <= 2, (dg.in_degree[node], dg.out_degree[node])
        assert dg.out_degree[node] <= 2, (dg.in_degree[node], dg.out_degree[node])

# これがうまくいったら、氷XVの部分秩序化に適用してみる。
