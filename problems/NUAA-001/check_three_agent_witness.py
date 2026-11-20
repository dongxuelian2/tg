"""Compare a Fano incidence graph's exact cut ratios with decentralized control LPs."""

from fractions import Fraction
from itertools import product
from random import Random

import numpy as np
from scipy.optimize import linprog


if __name__ == "__main__":
    # The three edge colors are the three cyclic matchings.
    edges = [(i, 7+(i+j) % 7) for i in range(7) for j in (0, 1, 3)]
    incidence = np.zeros((14, 21))
    for e, (i, j) in enumerate(edges):
        incidence[[i, j], e] = 1
    size = 1 << 14
    error = [Fraction(0)]*size
    first_bad = 15
    for mask in range(1, size):
        boundary = sum(((mask >> i) ^ (mask >> j)) & 1 for i, j in edges)
        worst = max(0, 1-Fraction(boundary, mask.bit_count()))
        # Dynamic program over submasks obtained by deleting one vertex.
        rest = mask
        while rest:
            bit = rest & -rest
            worst = max(worst, error[mask-bit])
            rest -= bit
        error[mask] = worst
        if worst:
            first_bad = min(first_bad, mask.bit_count())
    assert first_bad >= 6 and error[-1] == 1
    rng = Random(19)
    masks = [size-1]+[rng.randrange(1, size) for _ in range(149)]
    for mask in masks:
        selected = [v for v in range(14) if mask & (1 << v)]
        A = incidence[selected]
        goals = np.where(np.array(selected) < 7, 1., -1.)
        rows = np.c_[np.r_[A, -A], -np.ones(2*len(selected))]
        lp = linprog([0.]*21+[1.], A_ub=rows, b_ub=np.r_[goals, -goals],
                     bounds=[(-1, 1)]*21+[(0, None)])
        assert lp.success and abs(lp.fun-float(error[mask])) < 1e-9
    print("Minimum number of scenarios with positive error:", first_bad)
    print("16383 subset cut values computed; 150 independent LPs agreed")
    signs = np.array(list(product((-1, 0, 1), repeat=5)))
    signs = signs[np.any(signs, axis=1)]
    for _ in range(50):
        ed = [rng.sample(range(5), 2) for _ in range(8)]
        A = np.zeros((5, len(ed)))
        for j, ends in enumerate(ed):
            A[ends, j] = 1
        k = np.array([rng.uniform(.1, 2) for _ in ed])
        g = np.array([rng.uniform(-3, 3) for _ in range(5)])
        certificate = (signs@g-abs(signs@A)@k)/abs(signs).sum(axis=1)
        lp = linprog([0.]*len(ed)+[1.], A_ub=np.c_[np.r_[A,-A], -np.ones(10)],
                     b_ub=np.r_[g,-g], bounds=[(-x, x) for x in k]+[(0, None)])
        assert lp.success and abs(lp.fun-max(0, max(certificate))) < 1e-8
    print("50 general signed-subset comparisons agreed")
