"""Check guard elimination and asynchronous clauses by exhaustive policy replay."""

from itertools import product, combinations
from random import Random


def explicit(w, chosen, loss):
    seq, ok, baseline = w
    for count in range(min(len(seq), loss)+1):
        for erased in combinations(range(len(seq)), count):
            first = next((j for j, key in enumerate(seq)
                          if key in chosen and j not in erased), None)
            if not (baseline if first is None else ok[first]):
                return False
    return True


def prune(worlds, loss):
    keep = set(range(4))
    while True:
        remove = set()
        for seq, ok, _ in worlds:
            earlier = 0
            for a, value in zip(seq, ok):
                if a in keep:
                    if not value and earlier <= loss:
                        remove.add(a)
                    earlier += 1
        if not remove:
            break
        keep -= remove
    if all(base or sum(a in keep for a in seq) > loss for seq, _, base in worlds):
        return keep
    return None


def asynchronous(rng):
    # Both the next local observation and safety depend on prior mode changes.
    observations = {(s, i, active): rng.randrange(3) for s in range(3)
                    for i in range(2) for active in product((False, True), repeat=2)}
    safe = {(s, active): rng.random() < .82 for s in range(3)
            for active in product((False, True), repeat=2)}

    def replay(table=None, forced=None):
        active, earlier, switches = [False, False], set(), set()
        result = True
        for s in range(3):
            before = tuple(active)
            for i in range(2):
                if active[i]:
                    continue
                key = 3*i + observations[s, i, before]
                take = table[key] if table is not None else forced[i] == s
                if take:
                    switches.add(key)
                    active[i] = True
                else:
                    earlier.add(key)
            result &= safe[s, tuple(active)]
        return result, earlier, switches

    bad = []
    for forced in product(range(4), repeat=2):
        good, earlier, switches = replay(forced=forced)
        if not good:
            bad.append((earlier, switches))
    for table in product((False, True), repeat=6):
        sat = all(any(table[a] for a in earlier) or any(not table[a] for a in switches)
                  for earlier, switches in bad)
        assert sat == replay(table=table)[0]


if __name__ == "__main__":
    rng = Random(125)
    for _ in range(160):
        worlds = [([rng.randrange(4) for _ in range(5)],
                   [rng.random() < .75 for _ in range(5)], rng.random() < .5)
                  for _ in range(3)]
        loss = rng.randrange(3)
        wins = []
        for bits in product((False, True), repeat=4):
            chosen = {i for i, x in enumerate(bits) if x}
            if all(explicit(w, chosen, loss) for w in worlds):
                wins.append(chosen)
        greatest = prune(worlds, loss)
        assert (greatest is None) == (not wins)
        if wins:
            assert greatest in wins and greatest == set.union(*wins)
    for _ in range(100):
        asynchronous(rng)
    print("160 erasure games and 6400 asynchronous table replays matched")
