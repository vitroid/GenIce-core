import random

import numpy as np


def fourth_vector(v1, v2, v3):
    """
    Calculate the fourth vector from three tetrahedral vectors.
    """
    return -(v1 + v2 + v3)


# from genice
def quat2rotmat(q):
    a, b, c, d = q
    sp11 = a * a + b * b - (c * c + d * d)
    sp12 = -2.0 * (a * d + b * c)
    sp13 = 2.0 * (b * d - a * c)
    sp21 = 2.0 * (a * d - b * c)
    sp22 = a * a + c * c - (b * b + d * d)
    sp23 = -2.0 * (a * b + c * d)
    sp31 = 2.0 * (a * c + b * d)
    sp32 = 2.0 * (a * b - c * d)
    sp33 = a * a + d * d - (b * b + c * c)
    return np.array(
        [[sp11, sp12, sp13], [sp21, sp22, sp23], [sp31, sp32, sp33]]
    ).T


def compensate(vectors):
    """
    Arrange missing vectors for a tetrahedral node.
    """
    n = len(vectors)
    assert 0 < n <= 4, "At least one vector is required."

    if len(vectors) == 1:
        v0 = vectors[0]
        e0 = v0 / np.linalg.norm(v0)
        # ru is a random unit vector
        ru = np.random.random(3)
        ru /= np.linalg.norm(ru)
        # r2 is a random vector perpendicular to v0
        # It is the pivot to rotate v0
        r2 = np.cross(e0, ru)
        r2 /= np.linalg.norm(r2)
        # rotate v0 to make the second vector
        phih = np.radians(109.5 / 2)
        # quaternion for rotation
        a = np.cos(phih)
        b = -np.sin(phih) * r2[0]
        c = np.sin(phih) * r2[1]
        d = -np.sin(phih) * r2[2]
        # quaternion to rotation matrix
        R = quat2rotmat([a, b, c, d])
        # print(r2, r2@R)
        v1 = v0 @ R
        # VERIFY
        # e1 = v1 / np.linalg.norm(v1)
        vectors.append(v1)
    if len(vectors) == 2:
        v1, v2 = vectors
        y = v2 - v1
        z = v1 + v2
        ey = y / np.linalg.norm(y)
        ez = z / np.linalg.norm(z)
        ex = np.cross(ey, ez)
        L1 = np.linalg.norm(v1)
        L2 = np.linalg.norm(v2)
        L = (L1 + L2) / 2
        theta = np.radians(109.5 / 2)
        v3 = ex * L * np.sin(theta) - ez * L * np.cos(theta)
        # print(v1,v2,v3)
        # VERIFY
        # e1 = v1 / np.linalg.norm(v1)
        # e2 = v2 / np.linalg.norm(v2)
        # e3 = v3 / np.linalg.norm(v3)
        # Lx = np.linalg.norm(ex)
        # print(e1@e2, e2@e3, e3@e1, Lx)
        vectors.append(v3)
    if len(vectors) == 3:
        v4 = fourth_vector(*vectors)
        vectors.append(v4)
    # returns the newly added vectors
    return vectors[n:]


#################################################
# modified from GenIce
# https://github.com/vitroid/GenIce


def tip4p():
    """Interaction sites of TIP4P (?)
    Order: OHHM
    """
    L1 = 0.9572 / 10
    L2 = 0.15 / 10
    theta = np.radians(104.52)

    hy = L1 * np.sin(theta / 2)
    hz = L1 * np.cos(theta / 2)
    mz = L2
    sites = np.array(
        [[0.0, 0.0, 0.0], [0.0, hy, hz], [0.0, -hy, hz], [0.0, 0.0, mz]]
    )
    sites -= (sites[1] + sites[2]) / 18
    return sites, "OHHM"


def orient_water(vout1, vout2):
    """
    Prepare the rotation matrix from two outgoing vectors.
    """
    y = vout2 - vout1
    z = vout1 + vout2
    y /= np.linalg.norm(y)
    z /= np.linalg.norm(z)
    x = np.cross(y, z)
    return np.array([x, y, z])


def molecules_iter(dg, layout, watermodel=tip4p, max_iter=100, pbc=False):
    """
    Generate an atomic arrangements of water cluster fron the given digraph dg.

    A node of the graph can contain its positions as an attribute "pos".
    watermodel is a function that returns the intramolecular atomic coordinates.
    """
    for v in dg:
        # v is the label of a vertex
        nei_in = [x for x in dg.predecessors(v)]
        nei_out = [x for x in dg.successors(v)]
        # edge vectors from/to vertex v (origin is always v)
        v_in = [layout[i] - layout[v] for i in nei_in]
        v_out = [layout[i] - layout[v] for i in nei_out]
        if pbc:
            v_in = [v - np.floor(v + 0.5) for v in v_in]
            v_out = [v - np.floor(v + 0.5) for v in v_out]

        # Compensate missing tetrahedral vectors
        addv = compensate(v_in + v_out)
        miss_out = 2 - len(nei_out)
        v_out += random.sample(addv, miss_out)

        # Rotation matrix calculated from 2-in 2-out vectors.
        R = orient_water(*v_out)
        water, atomtypes = watermodel()
        yield water @ R + layout[v], atomtypes
