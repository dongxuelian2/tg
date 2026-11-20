"""Independent finite packet-history LP and exhaustive scheduling checks."""

from itertools import combinations, product
from random import Random

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


def tails(seq):
    return [sum(seq[i:]) for i in range(len(seq)+1)]


def bound(k, d, packets, reveal, target):
    K, D = tails(k), tails(d)
    return max([0]+[D[max([0]+[s for t, s in packets if t < j])]
                    + target*(j <= reveal) - K[j] for j in range(len(K))])


def lp(k, d, packets, reveal, target):
    indices, rows, cols, data, rhs = {}, [], [], [], []
    for wind in product((-1, 1), repeat=len(k)):
        drift = np.r_[0, np.cumsum(np.array(d)*wind)]
        for sign, side in product((-1, 1), repeat=2):
            row = len(rhs)
            rhs.append(side*(sign*target-drift[-1]))
            for j, gain in enumerate(k):
                reports = tuple(drift[s] for t, s in packets if t <= j)
                key = (j, reports, sign if j >= reveal else 0)
                col = indices.setdefault(key, len(indices)+1)
                rows.append(row)
                cols.append(col)
                data.append(side*gain)
            rows.append(row)
            cols.append(0)
            data.append(-1)
    arr = coo_matrix((data, (rows, cols)), shape=(len(rhs), len(indices)+1))
    result = linprog([1]+[0]*len(indices), A_ub=arr, b_ub=rhs,
                     bounds=[(0, None)]+[(-1, 1)]*len(indices))
    assert result.success
    return result.fun


def greedy(k, d, lag, delta):
    K, D = tails(k), tails(d)
    stamp, used = 0, 0
    while D[stamp] > delta:
        next_stamp = max((t-lag for t, cap in enumerate(K)
                          if cap >= D[stamp]-delta), default=-1)
        if next_stamp <= stamp:
            return None
        stamp, used = next_stamp, used+1
    return used


if __name__ == "__main__":
    rng = Random(29)
    for _ in range(64):
        n = 5
        k, d = ([rng.randint(0, 4) for _ in range(n)] for _ in range(2))
        packets = [(t, rng.randrange(t+1)) for t in range(1, n+1)
                   if rng.randrange(3)]
        reveal, target = rng.randrange(n+1), rng.randrange(5)
        assert abs(bound(k, d, packets, reveal, target)
                   - lp(k, d, packets, reveal, target)) < 1e-7
        lag, delta = rng.randrange(3), rng.randrange(12)
        available = [(t, t-lag) for t in range(lag+1, n+1)]
        minimum = next((m for m in range(len(available)+1)
                        if any(bound(k, d, p, 0, 0) <= delta
                               for p in combinations(available, m))), None)
        assert minimum == greedy(k, d, lag, delta)
    print("64 packet-history LPs and 64 exhaustive schedules agreed")
