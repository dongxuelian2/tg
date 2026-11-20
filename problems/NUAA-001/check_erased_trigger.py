"""Check the physical deadline formula against shared-action terminal LPs."""

from itertools import combinations_with_replacement, product
from random import Random

import numpy as np
from scipy.optimize import linprog


def paths(label, code, p, e, n):
    histories = [((), p, e)]
    for _ in range(n):
        after = []
        for h, tries, losses in histories:
            send = tries > 0 and code[label, h]
            after.append((h+(int(send),), tries-int(send), losses))
            if send and losses:
                after.append((h+(0,), tries-1, losses-1))
        histories = after
    return [h for h, _, _ in histories]


def compare(leaves, cap):
    result = max([0]+[abs(x)-sum(cap) for _, x in leaves])
    for (a, x), (b, y) in combinations_with_replacement(leaves, 2):
        first = next((j for j in range(len(cap)) if a[j] != b[j]), len(cap))
        result = max(result, abs(x-y)/2 - sum(cap[first:]))

    nodes = sorted({h[:j+1] for h, _ in leaves for j in range(len(cap))})
    rows, bound = [], []
    for h, x in leaves:
        uses = np.array([int(h[:len(node)] == node) for node in nodes])
        for sign in (-1, 1):
            rows.append([*list(sign*uses), -1])
            bound.append(sign*x)
    limits = [(-cap[len(h)-1], cap[len(h)-1]) for h in nodes] + [(0, None)]
    fit = linprog([0]*len(nodes)+[1], A_ub=rows, b_ub=bound, bounds=limits)
    assert fit.success and abs(result-fit.fun) < 1e-8
    return result


if __name__ == "__main__":
    rng = Random(711)
    for _ in range(96):
        n, p, e = 3, rng.randrange(4), rng.randrange(3)
        targets = [-2., .25, 1.5]
        keys = [(i, h) for i in range(len(targets)) for j in range(n)
                for h in product((0, 1), repeat=j)]
        code = {key: rng.randrange(2) for key in keys}
        leaves = [(h, x) for i, x in enumerate(targets)
                  for h in paths(i, code, p, e, n)]
        compare(leaves, [rng.uniform(0, 1) for _ in range(n)])
    keys = [(i, h) for i in range(2) for h in [(), (0,), (1,)]]
    cap = [.7, .4]
    for p, e in product(range(3), repeat=2):
        best = float("inf")
        for bits in product((0, 1), repeat=len(keys)):
            code = dict(zip(keys, bits))
            leaves = [(h, x) for i, x in enumerate((-1., 1.))
                      for h in paths(i, code, p, e, 2)]
            best = min(best, compare(leaves, cap))
        expect = max(0, 1-sum(cap[e:])) if p > e else 1
        assert abs(best-expect) < 1e-8
    print("96 random trees and 576 feedback-encoder LPs agreed")
