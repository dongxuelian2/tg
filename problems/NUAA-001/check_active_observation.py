"""Compare two safety games by independently enumerating physical observations."""

from functools import lru_cache
from random import Random


def verify(seed):
    rng = Random(seed)
    positions = range(6)
    signs = (-1, 1)
    bias, noise, horizon = 2, 1, 3
    signal = [rng.randrange(-4, 5) for _ in positions]
    safe = [{s for s in signs if rng.random() < .85} for _ in positions]
    targets = [{s for s in signs if rng.random() < .55} for _ in positions]
    # Every action has two possible physical successors, both publicly observed.
    graph = [[tuple(rng.choices(positions, k=2)) for _ in range(2)]
             for p in positions]

    @lru_cache(None)
    def physical(p, time, worlds):
        types = {s for s, _ in worlds}
        if not types <= safe[p]:
            return False
        observations = {}
        for s, b in worlds:
            for w in range(-noise, noise + 1):
                y = s*signal[p] + b + w
                observations.setdefault(y, set()).add((s, b))
        for branch in observations.values():
            if time == horizon:
                if not {s for s, _ in branch} <= targets[p]:
                    return False
            elif not any(all(physical(q, time + 1, frozenset(branch)) for q in action)
                         for action in graph[p]):
                return False
        return True

    @lru_cache(None)
    def reduced(p, time, low, high, known):
        types = (known,) if known else signs
        if not set(types) <= safe[p]:
            return False
        low, high = min(low, signal[p]), max(high, signal[p])
        ambiguous = high-low <= 2*noise and max(high, -low) <= bias+noise
        branches = (0,) if not known and ambiguous else types
        for branch in branches:
            if time == horizon:
                if not set((branch,) if branch else signs) <= targets[p]:
                    return False
            elif not any(all(reduced(q, time + 1, low, high, branch) for q in action)
                         for action in graph[p]):
                return False
        return True

    worlds = frozenset((s, b) for s in signs for b in range(-bias, bias + 1))
    a = physical(0, 0, worlds)
    b = reduced(0, 0, signal[0], signal[0], 0)
    assert a == b, seed
    return a


if __name__ == "__main__":
    wins = sum(verify(seed) for seed in range(500))
    print(f"500 games agreed: {wins} winning, {500-wins} losing")
