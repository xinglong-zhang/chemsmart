import numpy as np


def is_collinear(coords, tol=1e-5):
    """Check if three points are collinear using cross product method."""
    a = coords[0]
    b = coords[1]
    c = coords[2]
    vec1 = np.array(b) - np.array(a)
    vec2 = np.array(c) - np.array(a)
    cross_product = np.cross(vec1, vec2)
    return (
        np.linalg.norm(cross_product) < tol
    )  # If cross product is close to zero, points are collinear
