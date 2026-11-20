"""Compare scan-placement LPs to exhaustive placement and packet outputs."""

from itertools import product
from random import Random

import numpy as np
from scipy.optimize import linprog


def check_graph(rng):
    layers = [[0], [1, 2], [3, 4], [5, 6], [7]]
    edges = [(a, b) for lhs, rhs in zip(layers, layers[1:])
             for a in lhs for b in rhs if rng.random() < .8]
    paths = [p for p in product(*layers) if all(e in edges for e in zip(p, p[1:]))]
    costs = [rng.randrange(1, 7) for _ in range(6)]
    q = rng.randrange(1, 4)
    optimum = min((sum(c*b for c, b in zip(costs, mask))
                   for mask in product((0, 1), repeat=6)
                   if all(sum(mask[i-1] for i in p[1:-1]) >= q for p in paths)),
                  default=None)
    # In/out potentials for eight nodes; six selection columns at the end.
    rows = []
    for i in range(8):
        row = np.zeros(22)
        row[2*i+1], row[2*i] = 1, -1
        if 1 <= i <= 6:
            row[15+i] = -1
        rows.append(row)
    for i, j in edges:
        row = np.zeros(22)
        row[2*j], row[2*i+1] = 1, -1
        rows.append(row)
    bounds = [(0, q)]*16 + [(0, 1)]*6
    bounds[0], bounds[15] = (0, 0), (q, q)
    sol = linprog([0]*16+costs, A_ub=rows, b_ub=np.zeros(len(rows)), bounds=bounds)
    assert (optimum is None and not sol.success) or (
        sol.success and abs(sol.fun-optimum) < 1e-8)


def outputs(word, budget):
    return {v for v in product(*[(0, x) if x else (0,) for x in word])
            if sum(x != y for x, y in zip(word, v)) <= budget}


if __name__ == "__main__":
    rng = Random(711)
    for _ in range(120):
        check_graph(rng)
    words = list(product(range(3), repeat=3))
    for e in range(3):
        sets = {w:outputs(w, e) for w in words}
        for a, b in product(words, repeat=2):
            ab = sum(x != 0 and x != y for x, y in zip(a, b))
            ba = sum(y != 0 and x != y for x, y in zip(a, b))
            assert bool(sets[a] & sets[b]) == (max(ab, ba) <= e)
    print("120 placement LPs and 2187 transcript comparisons passed")
