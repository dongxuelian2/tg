"""Compare the exact wind bridge to two independent fixed-moment LPs per case."""

import numpy as np
from scipy.optimize import linprog


def upper(h, b, pos, vel, t):
    n = (h-vel/b)/2
    if n == 0:
        return b*t*t/2, []
    a = np.clip(h-n/2+(pos-b*h*h/2)/(2*b*n), 0, h-n)
    return b*t*t/2-b*(max(t-a, 0)**2-max(t-a-n, 0)**2), [a, a+n]


if __name__ == "__main__":
    rng = np.random.default_rng(312)
    for _ in range(120):
        h, b = rng.uniform(.3, 3, 2)
        vel = b*h*rng.uniform(-.98, .98)
        pos = h*vel/2+(b*b*h*h-vel*vel)/(4*b)*rng.uniform(-1, 1)
        t = h*rng.random()
        hi, high_cuts = upper(h, b, pos, vel, t)
        negative_low, low_cuts = upper(h, b, -pos, -vel, t)
        grid = np.unique([0, h, t, *high_cuts, *low_cuts])
        width = np.diff(grid)
        eq = np.vstack([width, width*(h-(grid[1:]+grid[:-1])/2)])
        cost = -np.diff(np.maximum(t-grid, 0)**2)/2
        for sign, predicted in [(1, -negative_low), (-1, -hi)]:
            lp = linprog(sign*cost, A_eq=eq, b_eq=[vel, pos], bounds=(-b, b))
            assert lp.success and abs(lp.fun-predicted) < 1e-7
    print("240 intermediate-position LP extrema agreed")
