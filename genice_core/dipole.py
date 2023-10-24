"""
Optimizes the orientations of directed paths to reduce the net dipole moment.
"""
from logging import getLogger, DEBUG
from typing import Union

import numpy as np
import networkx as nx


def net_polarization(
    dg: nx.DiGraph, vertexPositions: np.ndarray, isPeriodicBoundary: bool = False
) -> np.ndarray:
    """Net polarization (actually a vector sum) of a digraph

    Args:
        dg (nx.DiGraph): The digraph.
        vertexPositions (np.ndarray): Positions of the vertices.
        isPeriodicBoundary (bool, optional): If true, the vertex positions must be in fractional coordinate. Defaults to False.

    Returns:
        np.ndarray: net polarization
    """
    pol = np.zeros(3)
    for i, j in dg.edges():
        d = vertexPositions[j] - vertexPositions[i]
        if isPeriodicBoundary:
            d -= np.floor(d + 0.5)
        pol += d
    return pol


def minimize_net_dipole(
    paths: list[list],
    vertexPositions: np.ndarray,
    maxiter: int = 2000,
    isPeriodicBoundary: bool = False,
    targetPol: Union[np.ndarray, None] = None,
) -> list[list]:
    """Minimize the net polarization by flipping several paths.

    It is assumed that every vector has an identical dipole moment.

    Args:
        paths (list of list): List of directed paths. A path is a list of integer. A path with identical labels at first and last items are considered to be cyclic.
        pos (nx.ndarray[*,3]): Positions of the nodes.
        maxiter (int, optional): Number of random orientations for the paths. Defaults to 1000.
        pbc (bool, optional): If `True`, the positions of the nodes must be in the fractional coordinate system.
        target (np.ndarray, optional): Target value for the dipole-moment optimization.
    Returns:
        list of paths: Optimized paths.
    """
    logger = getLogger()

    if maxiter < 1:
        return paths

    if targetPol is None:
        targetPol = np.zeros_like(vertexPositions[0])

    # polarized chains and cycles. Small cycle of dipoles are eliminated.
    polarized = []

    dipoles = []
    for i, path in enumerate(paths):
        if isPeriodicBoundary:
            # vectors between adjacent vertices.
            displace = vertexPositions[path[1:]] - vertexPositions[path[:-1]]
            # PBC wrap
            displace -= np.floor(displace + 0.5)
            # total dipole along the chain (or a cycle)
            chain_pol = np.sum(displace, axis=0)
            # if it is large enough, i.e. if it is a spanning cycle,
            if chain_pol @ chain_pol > 1e-6:
                logger.debug(path)
                dipoles.append(chain_pol)
                polarized.append(i)
        else:
            # dipole moment of a path; NOTE: No PBC.
            if path[0] != path[-1]:
                # If no PBC, a chain pol is simply an end-to-end pol.
                chain_pol = vertexPositions[path[-1]] - vertexPositions[path[0]]
                dipoles.append(chain_pol)
                polarized.append(i)
    dipoles = np.array(dipoles)
    logger.debug(dipoles)

    pol_optimal = np.sum(dipoles, axis=0) - targetPol
    logger.info(f"init {np.linalg.norm(pol_optimal)} dipole")
    parity_optimal = np.ones(len(dipoles))
    for loop in range(maxiter):
        parity = np.random.randint(2, size=len(dipoles)) * 2 - 1
        net_pol = parity @ dipoles - targetPol
        if net_pol @ net_pol < pol_optimal @ pol_optimal:
            pol_optimal = net_pol
            parity_optimal = parity
            logger.info(f"{loop} {pol_optimal} dipole")
            if pol_optimal @ pol_optimal < 1e-10:
                logger.debug("Optimized.")
                break

    for i, dir in zip(polarized, parity_optimal):
        if dir < 0:
            paths[i] = paths[i][::-1]

    if logger.isEnabledFor(DEBUG):
        # assert the chains are properly inversed.
        dipoles = []
        for i, path in enumerate(paths):
            # dipole moment of a path; NOTE: No PBC.
            if path[0] != path[-1]:
                # If no PBC, a chain pol is simply an end-to-end pol.
                chain_pol = vertexPositions[path[-1]] - vertexPositions[path[0]]
                dipoles.append(chain_pol)
        dipoles = np.array(dipoles)

        pol = np.sum(dipoles, axis=0) - targetPol
        pol -= pol_optimal
        assert pol @ pol < 1e-20

    return paths
