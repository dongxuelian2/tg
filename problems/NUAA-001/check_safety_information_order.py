"""Finite independent checks of the information order and its physical witnesses."""

from itertools import combinations, product
from math import expm1

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linprog


def trace_constraints(reports):
    # Time step and command radius are one. Omit fixed zero initial commands.
    worlds, horizon = len(reports), len(reports[0])
    size = worlds * horizon
    eye = np.eye(size).reshape(worlds, horizon, size)
    rows, bounds = [], []

    def bound(vector, limit):
        rows.extend((vector, -vector))
        bounds.extend((limit, limit))

    for world in range(worlds):
        for t in range(horizon):
            previous = eye[world, t - 1] if t else 0
            bound(eye[world, t] - previous, 1)
    for left, right in combinations(range(worlds), 2):
        distance = 0
        for t in range(horizon):
            distance = max(distance, abs(reports[left][t] - reports[right][t]))
            bound(eye[left, t] - eye[right, t], distance)
    return np.array(rows), np.array(bounds), eye


def check_distance():
    fine = [[0, 0, 1, 5], [0, 2, 3, 4], [1, 1, 2, 2]]
    coarse = [[0, 0, 1, 1], [0, 0, 0, 1], [1, 0, 0, 2]]
    hulls = [trace_constraints(reports) for reports in (fine, coarse)]
    size = hulls[0][0].shape[1]
    largest_deficit = largest_error = 0
    for left, right in combinations(range(len(fine)), 2):
        distances = [0, 0]
        for t in range(len(fine[0])):
            for j, (reports, (rows, bounds, eye)) in enumerate(zip((fine, coarse), hulls)):
                prefix_gap = max(abs(a - b) for a, b in
                                 zip(reports[left][:t + 1], reports[right][:t + 1]))
                distances[j] = min(distances[j] + 2, prefix_gap)
                objective = eye[right, t] - eye[left, t]
                result = linprog(objective, A_ub=rows, b_ub=bounds,
                                 bounds=(None, None))
                assert result.success and abs(result.fun + distances[j]) < 1e-8
                if j == 0:
                    candidate = result.x
            deficit = max(0, distances[0] - distances[1])
            largest_deficit = max(largest_deficit, deficit)
            rows, bounds, _ = hulls[1]
            identity = np.eye(size)
            constraints = np.vstack((np.c_[rows, np.zeros(len(rows))],
                                     np.c_[identity, -np.ones(size)],
                                     np.c_[-identity, -np.ones(size)]))
            rhs = np.r_[bounds, candidate, -candidate]
            projection = linprog(np.r_[np.zeros(size), 1], A_ub=constraints,
                                 b_ub=rhs, bounds=(None, None))
            assert projection.success
            largest_error = max(largest_error, projection.fun)
    assert largest_deficit > 0 and abs(largest_error - largest_deficit / 2) < 1e-8
    return largest_deficit, largest_error


def check_motor():
    def dynamics(time, state):
        thrust_increment, speed, _ = state
        return 10 * (4 - thrust_increment), thrust_increment, speed

    result = solve_ivp(dynamics, (0, 0.2), (0, 0, 0), rtol=1e-11, atol=1e-13)
    assert result.success
    radius = result.y[2, -1]
    assert abs(radius + 0.04 * expm1(-2)) < 1e-12
    return radius


def check_protocols():
    for count in range(1, 6):
        worlds = range(2 * count)
        for size in range(2 * count + 1):
            for subset in combinations(worlds, size):
                available = any(not {2 * j, 2 * j + 1}.issubset(subset)
                                for j in range(count))
                assert available == (size < len(worlds))
    for n in (3, 5, 7):
        edges = [(j, (j + 1) % n) for j in range(n)]
        for missing in [None, *range(2 * n)]:
            active = [edge for j, edge in enumerate(edges)
                      if missing is None or missing // 2 != j]
            feasible = any(all(bits[a] != bits[b] for a, b in active)
                           for bits in product(range(2), repeat=n))
            assert feasible == (missing is not None)


if __name__ == "__main__":
    print("trace deficit and projection error:", check_distance())
    print("critical tube radius (m):", check_motor())
    check_protocols()
    print("finite protocol witnesses passed")
