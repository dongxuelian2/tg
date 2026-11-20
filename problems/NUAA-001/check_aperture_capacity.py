"""Independent integer checks of square covers and moving-wall intersections."""

from functools import lru_cache
from random import Random


def cover_number(points, side):
    masks = {sum(1 << i for i, (x, y) in enumerate(points)
                 if a <= x <= a+side and b <= y <= b+side)
             for a, _ in points for _, b in points}

    @lru_cache(None)
    def solve(rest):
        if not rest:
            return 0
        first = rest & -rest
        return 1 + min(solve(rest & ~mask) for mask in masks if mask & first)
    return solve((1 << len(points))-1)


def wall_number(points, side):
    # Input units are integers. Multiply coordinates by eight for exact walls.
    gaps = [abs(p[j]-q[j])-side for p in points for q in points for j in (0, 1)
            if abs(p[j]-q[j]) > side]
    gap = min([side, *gaps])
    halfwidth = 4*side+gap
    rectangles = [(8*x-halfwidth, 8*y-halfwidth, 8*x+halfwidth, 8*y+halfwidth)
                  for x, y in points]

    def place(i, classes):
        if i == len(points):
            return len(classes)
        r = rectangles[i]
        best = place(i+1, classes + [r])
        for j, c in enumerate(classes):
            joint = (max(r[0], c[0]), max(r[1], c[1]),
                     min(r[2], c[2]), min(r[3], c[3]))
            if joint[0] <= joint[2] and joint[1] <= joint[3]:
                best = min(best, place(i+1, classes[:j] + [joint] + classes[j+1:]))
        return best
    return place(0, [])


if __name__ == "__main__":
    pentagon = [(1, 0), (3, 0), (0, 1), (4, 1), (2, 3)]
    assert cover_number(pentagon, 2) == wall_number(pentagon, 2) == 3
    rng = Random(19)
    for _ in range(120):
        points = [(rng.randrange(-8, 9), rng.randrange(-8, 9)) for _ in range(7)]
        side = rng.randrange(1, 8)
        assert cover_number(points, side) == wall_number(points, side)
    print("Five-cycle and 120 point families agreed")
