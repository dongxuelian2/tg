"""Compare interval propagation with exhaustive lattice obstacle histories."""

from random import Random
from math import sqrt


def update(regions, travel, left, right):
    merged = []
    for a, b in regions:
        a, b = a-travel, b+travel
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    out = []
    for a, b in merged:
        if a <= left:
            out.append((a, min(b, left)))
        if b >= right:
            out.append((max(a, right), b))
    return out


if __name__ == "__main__":
    # Units of one tenth, to avoid floating-point tests at open scan boundaries.
    first = update([(-50, 50)], 0, -24, 14)
    second = update(first, 10, -2, 22)
    assert second == [(-60, -14), (22, 60)]
    assert [(a-10, b+10) for a, b in second] == [(-70, -4), (12, 70)]
    rng = Random(21)
    for _ in range(300):
        regions, feasible = [(-6, 6)], set(range(-6, 7))
        for _ in range(5):
            travel = rng.randrange(4)
            left, right = sorted(rng.sample(range(-10, 11), 2))
            feasible = {y for x in feasible for y in range(x-travel, x+travel+1)
                        if not left < y < right}
            regions = update(regions, travel, left, right)
            assert feasible == {x for a, b in regions for x in range(a, b+1)}
    required = min(max(sqrt((x+1)**2+1), sqrt((x-1)**2+1))
                   for x in [i/1000 for i in range(-1000, 1001)])
    assert abs(required-sqrt(2)) < 1e-12
    print("1500 scan updates, late-report example, and speed threshold agreed")
