"""Compare local correction LPs against graph potentials and nested-width bounds."""

import numpy as np
from scipy.optimize import linprog


def fit(rows, intervals, bounds):
    A = np.array(rows)
    ends = np.array(intervals)
    C = np.c_[np.r_[A, -A], -np.ones(len(rows)*2)]
    ans = linprog([0]*len(bounds)+[1], A_ub=C, b_ub=np.r_[ends[:, 0], -ends[:, 1]],
                  bounds=bounds+[(0, None)])
    assert ans.success
    return ans.fun


def private(rng):
    k = rng.uniform(.1, 2, 6)
    edges = [(0, i+1, x, 0) for i, x in enumerate(k)]
    edges += [(v, u, c, m) for u, v, c, m in edges]
    rows, ends = [], []
    for s in range(3):
        for t in range(3, 6):
            if rng.random() < .3:
                continue
            g = rng.normal()
            d = rng.uniform(0, .4)
            edges.extend([(t+1, s+1, g-d, 1), (s+1, t+1, -g-d, 1)])
            row = np.zeros(6)
            row[[s, t]] = 1
            rows.append(row)
            ends.append([g-d, g+d])
    if not rows:
        return
    optimum = fit(rows, ends, [(-x, x) for x in k])
    for delta in [optimum+1e-6, optimum-1e-6]:
        dist = np.zeros(7)
        for _ in range(6):
            for u, v, c, m in edges:
                dist[v] = min(dist[v], dist[u]+c+m*delta)
        feasible = all(dist[v] <= dist[u]+c+m*delta+1e-10 for u, v, c, m in edges)
        assert feasible == (delta > optimum)


def nested(rng):
    K = rng.uniform(.1, 1, 3)
    targets = rng.normal(size=(2, 2, 2, 2))
    exact = max(0., np.max(abs(targets))-sum(K))
    bounds = []
    for j, k in enumerate(K, start=1):
        bounds.extend([(-k, k)]*2**j)
        vals = targets.reshape(2**j, -1)
        exact = max(exact, max(np.ptp(vals, axis=1))/2-sum(K[j:]))
    rows = []
    for x in range(16):
        row = np.zeros(14)
        for j in range(1, 4):
            row[2**j-2+(x >> (4-j))] = 1
        rows.append(row)
    assert abs(fit(rows, [(x, x) for x in targets.flat], bounds)-exact) < 1e-8


if __name__ == "__main__":
    rng = np.random.default_rng(317)
    for _ in range(150):
        private(rng)
        nested(rng)
    print("150 private graph comparisons and 150 nested LP checks passed")
