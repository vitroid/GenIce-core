from logging import DEBUG, INFO, basicConfig, getLogger

import networkx as nx

# import py3Dmol
from genice2.genice import GenIce
from genice2.plugin import Format, Lattice

from genice_core import ice_graph

logger = getLogger()
basicConfig(level=INFO)

lattice = Lattice("1h")
formatter = Format("raw", stage=(1, 2))
raw = GenIce(lattice, signature="Ice Ih", rep=(3, 3, 3)).generate_ice(formatter)

# graph is the topology of the hydrogen-bond network
g = nx.Graph(raw["graph"])
# reppositions contains the positions of CoM of water in fractional coordinate
layout = raw["reppositions"]
# repcell is the cell matrix (transposed)
cell = raw["repcell"]

fixed = nx.DiGraph(
    [
        [0, 4],
    ]
)

# set orientations of the hydrogen bonds.
dg = ice_graph(
    g,
    vertexPositions=layout,
    isPeriodicBoundary=True,
    dipoleOptimizationCycles=100,
    fixedEdges=fixed,
)
