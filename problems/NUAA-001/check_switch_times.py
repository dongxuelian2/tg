"""Independent exhaustive checks on small private timing and repeated handover games."""
from itertools import product
from math import inf
from random import Random


def stopping(t, u, first_report):
    # Evaluate the local stop/no-stop and action history at each tick.
    return all((t <= k) == (u <= k) for k in range(min(first_report, max(t, u)+1)))


def crossing_check():
    count = 0
    for g, ra, rb in product((1, 2), range(5), range(5)):
        answer = g + min(max(ra, rb), g + min(ra, rb))
        best = inf
        for a1, b1, a2, b2 in product(range(answer+1), repeat=4):
            if b1-a1 < g or a2-b2 < g:
                continue
            if stopping(a1, a2, ra) and stopping(b1, b2, rb):
                best = min(best, max(a1, a2, b1, b2))
        assert best == answer
        count += 1
    return count


def closure_check():
    rng = Random(20261119)
    for _ in range(120):
        n = 3
        lower = [rng.randrange(3) for _ in range(n)]
        upper = [rng.randrange(3,7) for _ in range(n)]
        delay = rng.randrange(5)
        edges = [(rng.randrange(n), rng.randrange(n), rng.randrange(-2,4), rng.randrange(2))
                 for _ in range(5)]

        def apply(t, value, capped):
            return min(t, max(0,value)+delay) if capped else t+value

        feasible = [t for t in product(*(range(a,b+1) for a,b in zip(lower,upper)))
                    if all(t[v] >= apply(t[u], c, cap) for u,v,c,cap in edges)]
        earliest = lower.copy()
        while all(a <= b for a,b in zip(earliest, upper)):
            prev = earliest.copy()
            for u,v,c,cap in edges:
                earliest[v] = max(earliest[v],apply(prev[u],c,cap))
            if earliest == prev:
                break
        assert bool(feasible) == all(a<=b for a,b in zip(earliest,upper))
        if feasible:
            assert tuple(earliest) == tuple(map(min,zip(*feasible)))


if __name__ == "__main__":
    print("crossing instances:", crossing_check())
    closure_check()
    print("random systems checked:", 120)
