"""Exhaust arbitrary integer encoder cells and controls on small two-stage games."""

from functools import lru_cache
from itertools import product


def brute(prior, stages, delta, budget):
    @lru_cache(None)
    def go(states, j, p):
        if j == len(stages):
            return all(abs(x) <= delta for x in states)
        k, d, m = stages[j]
        costs = [0] + ([1]*m if p else [])
        for labels in product(range(len(costs)), repeat=len(states)):
            for a, cost in enumerate(costs):
                cell = [x for x, label in zip(states, labels) if label == a]
                if cell and not any(go(tuple(sorted({x+q+w for x in cell
                                    for w in range(-d, d+1)})), j+1, p-cost)
                                    for q in range(-k, k+1)):
                    break
            else:
                return True
        return False
    return go(tuple(prior), 0, budget)


def lattice_table(stages, delta, budget):
    envelope = delta
    counts = [2*delta+1]*(budget+1)
    for k, d, m in reversed(stages):
        envelope += k-d
        after = [max(0, count-2*d) for count in counts]
        counts = [max(0, min(2*envelope+1, after[p]+(m*after[p-1] if p else 0)))
                  for p in range(budget+1)]
    return envelope, counts[budget]


if __name__ == "__main__":
    trials = 0
    for a, b in product([(1, 0, 1), (0, 1, 1), (1, 1, 1), (1, 0, 2)], repeat=2):
        for p in range(3):
            env, count = lattice_table([a, b], 1, p)
            for left, right in [(-3, -2), (-2, 1), (-1, 0), (0, 0), (1, 2)]:
                prior = range(left, right+1)
                expect = len(prior) <= count and max(map(abs, prior)) <= env
                assert brute(prior, [a, b], 1, p) == expect, (a, b, p, prior)
                trials += 1
    print(f"{trials} integer games matched unrestricted encoder enumeration")
