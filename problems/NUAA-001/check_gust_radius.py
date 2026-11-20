"""Solve finite causal games by explicit wind histories and independent linear programs."""

from itertools import product

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


def solve(capacity, gust, reveal, target):
    count = len(gust)
    keys, rows, cols, data, bound = {}, [], [], [], []
    for winds in product((-1, 1), repeat=count):
        drift = np.dot(winds, gust)
        for mode in (-1, 1):
            for side in (-1, 1):
                row = len(bound)
                bound.append(side*(mode*target-drift))
                for j, gain in enumerate(capacity):
                    key = (winds[:j+1], mode if j >= reveal else 0)
                    column = keys.setdefault(key, len(keys)+1)
                    rows.append(row)
                    cols.append(column)
                    data.append(side*gain)
                rows.append(row)
                cols.append(0)
                data.append(-1)
    constraints = coo_matrix((data, (rows, cols)), shape=(len(bound), len(keys)+1))
    cost = np.zeros(len(keys)+1)
    cost[0] = 1
    answer = linprog(cost, A_ub=constraints, b_ub=bound,
                     bounds=[(0, None)] + [(-1, 1)]*len(keys))
    assert answer.success
    return answer.fun


if __name__ == "__main__":
    trials = 0
    for c, b in [(1., .15), (1., 2.), (.8, 0.)]:
        s = np.linspace(0, 1, 6)
        capacity = c*np.diff(s**5)[::-1]/5
        gust = b*np.diff(s**2)[::-1]/2
        for reveal in range(6):
            for target in (0., .04, .2):
                E = np.maximum(gust-capacity, 0).sum()
                S = np.maximum(capacity[reveal:]-gust[reveal:], 0).sum()
                assert abs(solve(capacity, gust, reveal, target) -
                           (E+max(0., target-S))) < 1e-8
                trials += 1
    rng = np.random.default_rng(18)
    for _ in range(30):
        capacity, gust = rng.uniform(0, .5, (2, 5))
        reveal = rng.integers(0, 6)
        target = rng.uniform(0, .8)
        suffix = np.r_[np.cumsum((gust-capacity)[::-1])[::-1], 0]
        exact = max(max(suffix[reveal:]), target+max(suffix[:reveal+1]))
        assert abs(solve(capacity, gust, reveal, target) - exact) < 1e-8
        trials += 1
    print(f"{trials} causal wind-tree LPs matched")
