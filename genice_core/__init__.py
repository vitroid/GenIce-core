"""
.. include:: ../README.md
"""

"""
Optimizes the orientations of directed paths to reduce the net dipole moment.
"""
import numpy as np
import networkx as nx
from genice_core.topology import noodlize, simple_paths
from genice_core.dipole import minimize_net_dipole, vector_sum
from typing import Union


def ice_graph(
    g: nx.Graph,
    vertexPositions: Union[np.ndarray, None] = None,
    isPeriodicBoundary: bool = False,
    dipoleOptimizationCycles: int = 0,
    fixed: Union[nx.DiGraph, None] = nx.DiGraph(),
) -> nx.DiGraph:
    """Make a digraph that obeys the ice rules.

    A new algorithm suggested by Prof. Sakuma, Yamagata University.

    Args:
        g (nx.Graph): A ice-like undirected graph.
        vertexPositions (Union[nx.ndarray, None], optional): Positions of the vertices. Defaults to None.
        isPeriodicBoundary (bool, optional): If True, the positions are considered to be in the fractional coordinate system. Defaults to False.
        dipoleOptimizationCycles (int, optional): Number of iterations to reduce the net dipole moment. Defaults to 0 (no iteration).
        fixed (Union[nx.DiGraph, None], optional): A digraph made of edges whose directions are fixed. All edges in fixed must also be included in g. Defaults to an empty graph.

    Returns:
        nx.DiGraph: An ice graph.
    """
    # Divide the graph into noodle graph
    divg = noodlize(g, fixed)

    # Simplify paths ( paths with least crossings )
    paths = list(simple_paths(len(g), divg))

    # arrange the orientations here if you want to balance the polarization
    if vertexPositions is not None:
        if fixed is not None:
            # Set the targetPol in order to cancel the polarization in fixed.
            targetPol = -vector_sum(fixed, vertexPositions, isPeriodicBoundary)
        else:
            targetPol = np.zeros(3)

        paths = minimize_net_dipole(
            paths,
            vertexPositions,
            isPeriodicBoundary=isPeriodicBoundary,
            dipoleOptimizationCycles=dipoleOptimizationCycles,
            targetPol=targetPol,
        )

    # paths to digraph
    dg = nx.DiGraph([(i, j) for path in paths for i, j in zip(path, path[1:])])
    return dg
