"""Compare the history metric to a graph and the analog maneuver bound to an LP."""
from itertools import product
import numpy as np
from scipy.integrate import quad
from scipy.optimize import linprog
from scipy.sparse.csgraph import shortest_path


def history_check():
    horizon = 4
    nodes = [w for k in range(horizon+1) for w in product((-1, 0, 1), repeat=k)]
    index = {w: j for j, w in enumerate(nodes)}
    cost = np.full((len(nodes), len(nodes)), np.inf)
    np.fill_diagonal(cost, 0)
    for j, w in enumerate(nodes):
        if w:
            cost[j, index[w[:-1]]] = cost[index[w[:-1]], j] = 0.7
        for k, v in enumerate(nodes):
            if len(v) == len(w) and w:
                cost[j, k] = (0.6+len(w)/9)*max(abs(a-b) for a, b in zip(w, v))
    paths = shortest_path(cost, directed=False)
    worst = 0
    for i, w in enumerate(nodes):
        for j, v in enumerate(nodes):
            explicit = 0.7*(len(w)+len(v))
            for r in range(1, min(len(w), len(v))+1):
                dr = max(abs(w[k]-v[k]) for k in range(r))
                explicit = min(explicit, 0.7*(len(w)+len(v)-2*r)+(0.6+r/9)*dr)
            worst = max(worst, abs(explicit-paths[i, j]))
    print("prefix pairs, maximum distance error:", len(nodes)**2, worst)
    assert worst < 1e-12


def maneuver_check():
    # The signal first becomes distinguishable at a grid point.
    H, D, k, rho, mass = 1.0, 1.8, 0.5, 0.125, 1.3
    n = 40
    t = np.linspace(0, H, n+1)
    weights = np.zeros(n+1)
    for j in range(n):
        a, b = t[j:j+2]
        for idx, shape in [(j, lambda s: (b-s)/(b-a)), (j+1, lambda s: (s-a)/(b-a))]:
            weights[idx] += quad(lambda s: -np.expm1(-D*(H-s))*shape(s), a, b)[0]
    for A, L in [(0.7, 0.4), (0.15, 0.9)]:
        # Unknowns are the two separate cumulative commands and the common guaranteed value.
        size = 2*(n+1)+1
        rows, rhs = [], []
        for side in range(2):
            offset = side*(n+1)
            for j in range(n):
                row = np.zeros(size)
                row[offset+j:offset+j+2] = [-1, 1]
                rows.extend([row, -row])
                rhs.extend([A*H/n]*2)
            row = np.zeros(size)
            row[offset:offset+n+1] = (-1 if side == 0 else 1)*4*weights/mass
            row[-1] = 1
            rows.append(row)
            rhs.append(0)
        for j in range(n+1):
            row = np.zeros(size)
            row[j], row[n+1+j] = 1, -1
            rows.extend([row, -row])
            rhs.extend([2*L*max(0, k*t[j]-rho)]*2)
        bounds = [(None, None)]*size
        bounds[0] = bounds[n+1] = (0, 0)
        result = linprog(np.r_[np.zeros(size-1), -1], A_ub=rows, b_ub=rhs, bounds=bounds)
        T = H-rho/k
        expected = 4*min(A,L*k)/mass*(T*T/2-T/D-np.expm1(-D*T)/D**2)
        assert result.success and abs(-result.fun-expected) < 1e-9
        print("guaranteed displacement, LP / formula:", -result.fun, expected)


if __name__ == "__main__":
    history_check()
    maneuver_check()
