"""Exhaust all pairs of four-label gate stages with two scheduled binary packets."""

from itertools import combinations, product


edges = list(combinations(range(4), 2))
vertices = range(4)


def compatible_edges(labels):
    return sum(1 << k for k, (i, j) in enumerate(edges) if labels[i] != labels[j])


def direct_protocols():
    # At each gate a packet history must determine a single safe preposition.
    plans = set()
    for final in product(range(4), repeat=4):
        first = [word // 2 for word in final]
        plans.add((compatible_edges(first), compatible_edges(final)))
    return plans


def can_split(group, edge_mask):
    for bits in product(range(2), repeat=len(group)):
        paint = dict(zip(group, bits))
        if all(not (edge_mask >> j & 1) or
               i not in paint or k not in paint or paint[i] != paint[k]
               for j, (i, k) in enumerate(edges)):
            return True
    return False


def feasible(first_mask, last_mask):
    for colors in product(range(2), repeat=4):
        if first_mask & ~compatible_edges(colors):
            continue
        if all(can_split([i for i in vertices if colors[i] == color], last_mask)
               for color in (0, 1)):
            return True
    return False


if __name__ == "__main__":
    plans = direct_protocols()
    wins = 0
    for early, late in product(range(64), repeat=2):
        answer = feasible(early, late)
        assert answer == any(not (early & ~a or late & ~b) for a, b in plans)
        wins += answer
    assert not feasible(7, 56) and feasible(0, 63)
    print(f"4096 games agreed: {wins} winning, {4096-wins} losing")
