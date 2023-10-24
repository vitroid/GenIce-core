"""
固定辺の処理。

とりあえず試験。
"""
import pairlist as pl
import numpy as np
import networkx as nx
import yaplotlib as yap
from collections import defaultdict
import random
from logging import getLogger, basicConfig, INFO, DEBUG
import genice_core
from genice_core.dipole import net_polarization


def eineout_to_digraph(ein, eout):
    dg = nx.DiGraph()
    for i in eout:
        for j in eout[i]:
            assert i in ein[j]
            if i >= 0 and j >= 0:
                dg.add_edge(i, j)
    return dg


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
basicConfig(level=DEBUG)


random.seed(0)

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
nodes = list(g)
for i in nodes:
    room = 4 - len(g[i])
    for j in range(room):
        g.add_edge(i, -j - 1)

# 固定矢印の配列
ein = defaultdict(list)
eout = defaultdict(list)

logger.debug((pos[16], list(g[18])))
logger.debug((pos[27], list(g[29])))
# assert False
# Dope ions
perimeter = []
NH4 = 18
neis = list(g[NH4])
for nei in neis:
    eout[NH4].append(nei)
    ein[nei].append(NH4)
    perimeter.append(nei)
    g.remove_edge(NH4, nei)

F = 29
neis = list(g[F])
for nei in neis:
    ein[F].append(nei)
    eout[nei].append(F)
    perimeter.append(nei)
    g.remove_edge(F, nei)


# perimeterのノードは、さらに辺を固定する必要があるかもしれない。
while len(perimeter) > 0:
    draw(pos, g, eineout_to_digraph(ein, eout))
    node = perimeter.pop(0)
    nnei = len(g[node])
    logger.info((node, g[node], ein[node], eout[node]))
    assert len(g[node]) + len(ein[node]) + len(eout[node]) == 4
    # 状況はこれで把握できた。いくつ固定すべきか、どうやって判定するのか?
    # g上の結合数によらず、水分子は本来4本の水素結合を持っている。
    # einが1本固定されたら、eoutも同じ数だけ固定されねばならない。
    # ただし、その一本はdanglingかもしれない。
    # danglingの場合には-1を指定する。
    if len(ein[node]) == len(eout[node]):
        continue
    elif len(ein[node]) > len(eout[node]):
        neis = list(g[node])
        # logger.info(neis)
        e = random.choice(neis)
        eout[node].append(e)
        ein[e].append(node)
        g.remove_edge(node, e)
        if e >= 0:
            perimeter.append(e)
    else:  # len(ein[node]) < len(eout[node]):
        neis = list(g[node])
        e = random.choice(neis)
        ein[node].append(e)
        eout[e].append(node)
        g.remove_edge(node, e)
        if e >= 0:
            perimeter.append(e)


# これで残った無向グラフにgenice-coreして、ちゃんとice ruleが満足されることを示す。

# Fix部分の双極子モーメント和。
totalPol = np.zeros(3)
for i in eout:
    if i >= 0:
        for j in eout[i]:
            if j >= 0:
                d = pos[j] - pos[i]
                # PBC process here
                totalPol += d
logger.debug(f"{totalPol} net dipole before optimization")

# ダミーノードを除去する
try:
    g.remove_node(-1)
    g.remove_node(-2)
    g.remove_node(-3)
    g.remove_node(-4)
except:
    pass

# 固定されなかった部分を作る。固定部分の分極をキャンセルするように。
dg = genice_core.ice_graph(
    g, vertexPositions=pos, targetPol=-totalPol, dipoleOptimizationCycles=100
)
draw(pos, nx.Graph(), dg)


for i in eout:
    if i >= 0:
        for j in eout[i]:
            if j >= 0:
                dg.add_edge(i, j)


totalPol = net_polarization(dg, pos)
logger.debug(f"{totalPol} net dipole after optimization")


draw(pos, nx.Graph(), dg)


# できたっぽい。
# いや、できてない。
# 2fixで2本残ったgの頂点は「縦割り」してはいけないが、ただの2結合頂点は縦割りもOK。この違いをプログラムするのは難しいぞ。
#
