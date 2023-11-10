"""
.. include:: ../README.md
"""

"""
Optimizes the orientations of directed paths to reduce the net dipole moment.
"""
import numpy as np
import networkx as nx
from genice_core.topology import noodlize, split_into_simple_paths, balance
from genice_core.dipole import optimize, vector_sum
from typing import Union
from logging import getLogger


def ice_graph(
    g: nx.Graph,
    vertexPositions: Union[np.ndarray, None] = None,
    isPeriodicBoundary: bool = False,
    dipoleOptimizationCycles: int = 0,
    fixedEdges: Union[nx.DiGraph, None] = nx.DiGraph(),
    hook=None,
) -> nx.DiGraph:
    """Make a digraph that obeys the ice rules.

    A new algorithm based on the suggestion by Prof. Sakuma, Yamagata University.

    Args:
        g (nx.Graph): A ice-like undirected graph.
        vertexPositions (Union[nx.ndarray, None], optional): Positions of the vertices. Defaults to None.
        isPeriodicBoundary (bool, optional): If True, the positions are considered to be in the fractional coordinate system. Defaults to False.
        dipoleOptimizationCycles (int, optional): Number of iterations to reduce the net dipole moment. Defaults to 0 (no iteration).
        fixed (Union[nx.DiGraph, None], optional): A digraph made of edges whose directions are fixed. All edges in fixed must also be included in g. Defaults to an empty graph.

    Returns:
        nx.DiGraph: An ice graph (fixed part is excluded).
    """
    logger = getLogger()

    logger.debug(g)
    logger.debug(fixedEdges)

    if fixedEdges is not None:
        balance(fixedEdges, g, hook=hook)

    # Divide the graph into noodle graph
    dividedGraph = noodlize(g, fixedEdges)

    # Simplify paths ( paths with least crossings )
    paths = list(split_into_simple_paths(len(g), dividedGraph))

    # arrange the orientations here if you want to balance the polarization
    if vertexPositions is not None:
        if fixedEdges is not None:
            # Set the targetPol in order to cancel the polarization in fixed.
            targetPol = -vector_sum(fixedEdges, vertexPositions, isPeriodicBoundary)
        else:
            targetPol = np.zeros(3)

        paths = optimize(
            paths,
            vertexPositions,
            isPeriodicBoundary=isPeriodicBoundary,
            dipoleOptimizationCycles=dipoleOptimizationCycles,
            targetPol=targetPol,
        )

    # paths to digraph
    dg = nx.DiGraph()
    for path in paths:
        nx.add_path(dg, path)

    return dg
