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

    logger.debug(g)
    logger.debug(fixedEdges)

    # 派生した環
    free_cycles = []

    # 指定された辺だけでなく、均衡するように拡大したグラフ。
    if fixedEdges.size() > 0:
        extendedFixedEdges = None
        while extendedFixedEdges is None:
            extendedFixedEdges, free_cycles = balance(fixedEdges, g, hook=hook)

        if logger.isEnabledFor(DEBUG):
            n_spanning = 0
            for cycle in free_cycles:
                chainPol = _dipole_moment_pbc(cycle, vertexPositions)
                if chainPol @ chainPol > 1e-6:
                    n_spanning += 1
            logger.debug(f"Spanning cycles: {n_spanning}")

    else:
        extendedFixedEdges = nx.DiGraph()

    # 今のところ、free_cyclesは活用できていない。
    # バランスの過程で大量の環を見付けているので、それを忘れてまたnoodlizeするのはもったいない。
    # noodlizeから除外すればいいのだが、ちょっと書きかたがややこしくなるので悩ましい。

    # Divide the graph into noodle graph
    dividedGraph = noodlize(g, extendedFixedEdges)

    # Simplify paths ( paths with least crossings )
    paths = list(split_into_simple_paths(len(g), dividedGraph))

    # arrange the orientations here if you want to balance the polarization
    if vertexPositions is not None:
        if fixedEdges is not None:
            # Set the targetPol in order to cancel the polarization in fixed.
            targetPol = -vector_sum(
                extendedFixedEdges, vertexPositions, isPeriodicBoundary
            )
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

    dg = nx.DiGraph(extendedFixedEdges)
    for path in paths:
        nx.add_path(dg, path)

    return dg
