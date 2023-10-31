"""
.. include:: ../README.md
"""

"""
Optimizes the orientations of directed paths to reduce the net dipole moment.
"""
import numpy as np
import genice_core.networky as nx
from genice_core.topology import noodlize, split_into_simple_paths, balance
from genice_core.dipole import optimize, vector_sum
from typing import Union
from logging import getLogger
import time


def ice_graph(
    g: nx.Graph,
    vertexPositions: Union[np.ndarray, None] = None,
    isPeriodicBoundary: bool = False,
    dipoleOptimizationCycles: int = 0,
    fixedEdges: Union[nx.DiGraph, None] = nx.DiGraph(),
    hook=None,
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
        nx.DiGraph: An ice graph (fixed part is excluded).
    """
    logger = getLogger()

    logger.debug(g)
    logger.debug(fixedEdges)

    if fixedEdges is not None:
        balance(fixedEdges, g, hook=hook)

    # Divide the graph into noodle graph
    now = time.time()
    dividedGraph = noodlize(g, fixedEdges)
    logger.info(f"{time.time() - now} sec noodlize()")

    now = time.time()
    # Simplify paths ( paths with least crossings )
    paths = list(split_into_simple_paths(len(g), dividedGraph))
    logger.debug(f"{[path for path in paths]} paths ############################")
    logger.info(f"{time.time() - now} sec split_into_simple_paths()")

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
