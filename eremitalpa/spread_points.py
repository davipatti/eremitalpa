from typing import Iterable
import numpy as np


def repulsion(points, k):
    force = np.zeros_like(points)
    for i, a in enumerate(points):
        for j, b in enumerate(points):
            if i != j:
                dist = a - b
                abs_force = abs(1 / dist**3)
                if dist < 0:
                    force[i] -= abs_force
                else:
                    force[i] += abs_force
    return force * k


def attraction(current, target, k):
    return (target - current) * k


def total(current, target, k_repel, k_attract):
    return repulsion(current, k_repel) + attraction(current, target, k_attract)


def spread_points(
    points: Iterable[float],
    tol: float = 1e-4,
    maxiter: int = 100,
    repel: float = 0.5,
    attract: float = 0.1,
) -> np.array:
    """
    Spread out 1D points. Imagines points are attracted to their starting poisitions
    (with `attract` force constant) and are repelled from other points with a force
    proportional to the inverse of their cubed distance (multiplied by `repel`).
    Iteratively updates poisitions until either `maxiter` is reached or the sum of the
    squared differences between positions in successive iterations falls below `tol`.
    """
    points = np.array(points)
    iterations = 0
    diff = np.inf
    orig = points
    # prevent points overlapping precisely, which would lead to infinite forces
    points = points + np.random.uniform(-0.005, 0.005, len(points))
    while diff > tol and iterations < maxiter:
        updated = points + total(points, target=orig, k_repel=repel, k_attract=attract)
        diff = ((points - updated) ** 2).sum()
        points = updated
        iterations += 1
    return updated
