from logging import DEBUG, INFO, basicConfig, getLogger

import networkx as nx

# import py3Dmol
from genice2.genice import GenIce
from genice2.plugin import Format, Lattice

from genice_core.gromacs import render
from genice_core.topology import ice_graph
from genice_core.water import tip4p

logger = getLogger()
basicConfig(level=DEBUG)

lattice = Lattice("1h")
formatter = Format("raw", stage=(1, 2))
raw = GenIce(lattice, signature="Ice Ih", rep=(3, 3, 3)).generate_ice(
    formatter
)

# graph is the topology of the hydrogen-bond network
g = nx.Graph(raw["graph"])
# reppositions contains the positions of CoM of water in fractional coordinate
layout = raw["reppositions"]
# repcell is the cell matrix (transposed)
cell = raw["repcell"]

# set orientations of the hydrogen bonds.
dg = ice_graph(g, pos=layout, pbc=True)

# put water molecules
gro = render(
    dg,
    layout @ cell,
    watermodel=tip4p,
    cell=f"{cell[0,0]} {cell[1,1]} {cell[2,2]}",
    pbc=True,
)
with open(f"save.gro", "w") as f:
    f.write(gro)

# # show
# view = py3Dmol.view()
# view.addModel(gro, "gro")
# view.setStyle({"stick": {}})
# view.addUnitCell()
# view.zoomTo()
# view.show()
