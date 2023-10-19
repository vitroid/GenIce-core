"""
.. include:: ../README.md
"""

"""
Optimizes the orientations of directed paths to reduce the net dipole moment.
"""
import numpy as np
import networkx as nx
from genice_core.topology import make_digraph, noodlize
from typing import Union


def ice_graph(
    g: nx.Graph,
    vertexPositions: Union[np.ndarray, None] = None,
    isPeriodicBoundary: bool = False,
    dipoleOptimizationCycles: int = 0,
) -> nx.DiGraph:
    """Make a digraph that obeys the ice rules.

    A new algorithm suggested by Prof. Sakuma, Yamagata University.

    Args:
        g (nx.Graph): A ice-like undirected graph.
        vertexPositions (Union[nx.array, None], optional): Positions of the vertices. Defaults to None.
        isPeriodicBoundary (bool, optional): If True, the positions are considered to be in the fractional coordinate system. Defaults to False.
        dipoleOptimizationCycles (int, optional): Number of iterations to reduce the net dipole moment. Default is 0 (no iteration).
    Returns:
        nx.DiGraph: An ice graph.
    """
    divg = noodlize(g)
    dg = make_digraph(
        g,
        divg,
        pos=vertexPositions,
        pbc=isPeriodicBoundary,
        dipoleOptimizationCycles=dipoleOptimizationCycles,
    )
    return dg
