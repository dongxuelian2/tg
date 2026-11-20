"""Exact rational conflict-graph checks; each edge has two physical witnesses."""

from fractions import Fraction as F
from math import ceil, floor


def colorable(edges, count):
    used = {}

    def search():
        if len(used) == len(edges):
            return True
        def forbidden(v):
            return {used[w] for w in edges[v] if w in used}
        v = max(edges.keys() - used.keys(),
                key=lambda v: (len(forbidden(v)), len(edges[v])))
        choices = set(range(min(count, len(set(used.values())) + 1))) - forbidden(v)
        for c in choices:
            used[v] = c
            if search():
                return True
            del used[v]
        return False

    return search()


def verify(ratio):
    low, high = F(1), ratio
    a, eta = (high+low)/8, (high-low)/8
    n = floor(ratio) + 1
    if ratio.denominator == 1:
        pts = set(range(n))
    else:
        s = ceil(1/(n*(ratio/(n-1)-1)))
        long = 1 + F(1, s*n)
        pts = {j*long for j in range(s*n+1)} | {F(1+j) for j in range(s*n+1)}
    radius = F(max(pts), 2)
    pts = {p-radius for p in pts}
    edges = {p: set() for p in pts}
    colors = ceil(ratio) + 1
    for p in pts:
        for q in pts:
            if low <= q-p <= high:
                z = (p+q)/2
                for x, sign in [(p, -1), (q, 1)]:
                    beta = (x+z)/2
                    u = x-beta-sign*a
                    v = z-beta+sign*a
                    assert abs(beta) <= radius and max(abs(u), abs(v)) <= eta
                edges[p].add(q)
                edges[q].add(p)
                assert floor(p/low) % colors != floor(q/low) % colors
    assert not colorable(edges, colors-1)
    return colors


if __name__ == "__main__":
    for ratio in map(F, ["5/4", "3/2", "2", "9/4", "5/2", "3", "7/2", "9/2"]):
        print(ratio, verify(ratio))
