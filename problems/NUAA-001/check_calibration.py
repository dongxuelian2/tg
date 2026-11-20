"""Independent enumeration of threshold trees on a rational threshold grid."""

from functools import lru_cache


def can_classify(radius, noise, separation, rounds):
    # Coordinate unit is 1/8. Retain the upper-end flag exactly.
    scale = 8
    radius, noise, separation = (scale * x for x in (radius, noise, separation))

    @lru_cache(None)
    def search(left, right, upper_closed, remaining):
        width = right - left
        if width < 0 or (width == 0 and not upper_closed):
            return True
        if width < 2 * separation or (width == 2 * separation and not upper_closed):
            return True
        if not remaining:
            return False
        for threshold in range(left - noise, right + noise + 1):
            cut = threshold + noise
            low = search(left, min(right, cut), upper_closed and right < cut,
                         remaining - 1)
            if low and search(max(left, threshold - noise), right, upper_closed,
                              remaining - 1):
                return True
        return False

    return search(-radius, radius, True, rounds)


if __name__ == "__main__":
    checks = 0
    for rounds in range(4):
        for noise, separation in [(0, 1), (1, 1), (2, 1), (1, 2)]:
            critical = separation + (2**rounds - 1) * max(separation - noise, 0)
            for radius in (critical - 1, critical, critical + 1):
                observed = can_classify(radius, noise, separation, rounds)
                assert observed == (radius < critical), (radius, noise, separation, rounds)
                checks += 1
    print(f"{checks} threshold-tree comparisons passed")
