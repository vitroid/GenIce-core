"""
.. include:: ../README.md
"""

"""
Optimizes the orientations of directed paths to reduce the net dipole moment.
"""
import numpy as np
import networkx as nx
from genice_core.topology import noodlize, split_into_simple_paths, balance
from genice_core.dipole import optimize, vector_sum, _dipole_moment_pbc
from typing import Union
from logging import getLogger, DEBUG


def ice_graph(
    g: nx.Graph,
    vertexPositions: Union[np.ndarray, None] = None,
    isPeriodicBoundary: bool = False,
    dipoleOptimizationCycles: int = 0,
    fixedEdges: nx.DiGraph = nx.DiGraph(),
    hook=None,
) -> nx.DiGraph:
    """Make a digraph that obeys the ice rules.

    A new algorithm based on the suggestion by Prof. Sakuma, Yamagata University.

    Args:
        g (nx.Graph): A ice-like undirected graph.
        vertexPositions (Union[nx.ndarray, None], optional): Positions of the vertices. Defaults to None.
        isPeriodicBoundary (bool, optional): If True, the positions are considered to be in the fractional coordinate system. Defaults to False.
        dipoleOptimizationCycles (int, optional): Number of iterations to reduce the net dipole moment. Defaults to 0 (no iteration).
        fixed (nx.DiGraph, optional): A digraph made of edges whose directions are fixed. All edges in fixed must also be included in g. Defaults to an empty graph.

    Returns:
        nx.DiGraph: An ice graph. (CHANGED. BE CAREFUL この変更の影響をもっとも受けるのがGenIce)
    """
    logger = getLogger()

    # derived cycles in extending the fixed edges.
    derivedCycles = []

    if fixedEdges.size() > 0:
        # balance fixed edges
        extendedFixedEdges = None
        while extendedFixedEdges is None:
            # It returns Nones when it fails to balance.
            extendedFixedEdges, derivedCycles = balance(fixedEdges, g, hook=hook)
    else:
        extendedFixedEdges = nx.DiGraph()

    # all edges that are oriented in balance()
    fullFixedEdges = nx.DiGraph(extendedFixedEdges)
    for cycle in derivedCycles:
        nx.add_path(fullFixedEdges, cycle)

    # Divide the remaining (unfixed) part of the graph into a noodle graph
    dividedGraph = noodlize(g, fullFixedEdges)

    # Simplify paths ( paths with least crossings )
    paths = list(split_into_simple_paths(len(g), dividedGraph)) + derivedCycles

    # arrange the orientations here if you want to balance the polarization
    if vertexPositions is not None:
        # Set the targetPol in order to cancel the polarization in the fixed part.
        targetPol = -vector_sum(extendedFixedEdges, vertexPositions, isPeriodicBoundary)

        paths = optimize(
            paths,
            vertexPositions,
            isPeriodicBoundary=isPeriodicBoundary,
            dipoleOptimizationCycles=dipoleOptimizationCycles,
            targetPol=targetPol,
        )

    # Combine everything together
    dg = nx.DiGraph(extendedFixedEdges)
    for path in paths:
        nx.add_path(dg, path)

    # # この時点で、pathsを検査しておく。
    # if logger.isEnabledFor(DEBUG):
    #     gg = nx.Graph(extendedFixedEdges)
    #     for path in paths:
    #         nx.add_path(gg, path)
    #     logger.debug(f"Size g {g.number_of_nodes()} {g.number_of_edges()}")
    #     logger.debug(f"Size gg {gg.number_of_nodes()} {gg.number_of_edges()}")
    #     assert g.number_of_edges() == gg.number_of_edges()
    #     e1 = set([(min(i, j), max(i, j)) for i, j in g.edges()])
    #     e2 = set([(min(i, j), max(i, j)) for i, j in gg.edges()])
    #     logger.debug(
    #         f"{sorted(list(e1 - e2))} edges only in original undirected graph."
    #     )
    #     logger.debug(f"{sorted(list(e2 - e1))} edges only in derived directed graph.")

    return dg
