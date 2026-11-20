"""Independent numerical falsification checks; run with Python, NumPy and SciPy."""

from itertools import combinations, product
from math import cos, sin, sqrt

import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.optimize import brentq, linprog, minimize


def one_switch(h, B, a=0.0, rate=0.0):
    # Integration variable y maps the suffix to [0,1], including q=h.
    def residual(q):
        def integrand(y):
            t = q + (h - q) * y
            angle = a + rate * t + B * (q * t - q**2 / 2)
            return y * (1 - y) * cos(angle)

        gain = B * (h - q) ** 2 * quad(integrand, 0, 1)[0]
        return gain - sin(a + rate * q + B * q**2 / 2) / 2

    return brentq(residual, 0, h, xtol=1e-15)


def support(h, B, a=0.0, rate=0.0):
    q = one_switch(h, B, a, rate)

    def state(t):
        moment = min(t, q)
        angle = a + rate * t + B * (moment * t - moment**2 / 2)
        return angle, 0.5 if t < q else 1.0

    def integrand(t, rate_gradient):
        angle, thrust = state(t)
        response = B * t * cos(angle) if rate_gradient else sin(angle)
        return (h - t) * thrust * response

    value, K = [quad(integrand, 0, h, args=(grad,), points=[q])[0]
                for grad in (False, True)]
    return value, K, q


def rotor_search(h, B, a, rate, rng, parts=32):
    # Optimize the two banks on each cell directly. No switch or sign restriction.
    nodes, weights = leggauss(6)
    step = h / parts
    times = (np.arange(parts)[:, None] + (1 + nodes) / 2) * step
    times = times.ravel()
    starts = np.arange(parts) * step
    age = np.maximum(times[:, None] - starts, 0)
    kernel = B * (age**2 - np.maximum(age - step, 0)**2) / 2
    quadrature = np.tile(weights * step / 2, parts) * (h - times)
    cells = np.repeat(np.arange(parts), len(nodes))

    def objective(banks):
        left, right = banks.reshape(2, parts)
        thrust = ((left + right) / 2)[cells]
        angles = a + rate * times + kernel @ (right - left)
        direct = np.bincount(cells, weights=quadrature * np.sin(angles)) / 2
        rotation = kernel.T @ (quadrature * thrust * np.cos(angles))
        return (-np.dot(quadrature * thrust, np.sin(angles)),
                -np.r_[direct - rotation, direct + rotation])

    values = []
    for _ in range(3):
        result = minimize(objective, rng.uniform(size=2 * parts), jac=True,
                          bounds=[(0, 1)] * (2 * parts), method="L-BFGS-B",
                          options={"gtol": 1e-13, "ftol": 1e-15, "maxiter": 1000})
        assert result.success, result.message
        values.append(-result.fun)
    return max(values)


def check_rotors():
    rng = np.random.default_rng(421)
    for B in (0.1, 1.0, 2 * np.pi / 3):
        value, K, q = support(1.0, B)
        # One-sided derivative stays in the specified short-horizon domain.
        eps = 1e-5
        derivative = (3 * value - 4 * support(1-eps, B)[0]
                      + support(1-2*eps, B)[0]) / (2*eps)
        assert abs(derivative - K) < 1e-8
        grid_value = rotor_search(1, B, 0, 0, rng)
        assert -1e-10 <= value - grid_value < 0.002 * value
        print(f"B={B:.6g} q={q:.9f} X={value:.9f} K={K:.9f} grid={grid_value:.9f}")
    # Independent small-angle limit: sin(z)~z, cos(z)~1 in the root equation.
    expected = sqrt(2) / (sqrt(2) + sqrt(3))
    assert abs(one_switch(1, 1e-4) - expected) < 1e-9
    # Check inherited pitch/rate without resetting them. The optimizer sees only
    # actual bank amplitudes, including potential counterrotation.
    B, a, rate, last = 1.9, 0.0, 0.0, 0.0
    for end in (1/16, 1/4, 1):
        h = end - last
        value, _, q = support(h, B, a, rate)
        alternate = rotor_search(h, B, a, rate, rng)
        assert -1e-11 <= value - alternate < value * 0.003
        a += rate * h + B * (q*h - q*q/2)
        rate += B*q
        last = end


def history_lp(agents):
    """Minimax terminal error with one shared unknown sign for all agents.

    Tuple fields: cells, lag in cells, wind coefficient, control coefficient,
    signed formation weight, blind interval duration. No cancellation law is
    supplied to the solver. Each local policy has one free real per history.
    """
    histories = [list(product((-1, 1), repeat=n)) for n, *_ in agents]
    columns = {}
    constraints, constants = [], []
    for worlds in product(*histories):
        row, wind = {}, 0
        for i, ((n, lag, p, k, weight, duration), world) in enumerate(zip(agents, worlds)):
            dt = duration / n
            wind += weight * dt * p * sum(world)
            for j in range(n):
                observed = world[:max(0, j + 1 - lag)]
                key = i, j, observed
                col = columns.setdefault(key, len(columns))
                row[col] = weight * dt * k
        constraints.append(row)
        constants.append(wind)
    matrix = np.zeros((len(constraints), len(columns)))
    for i, row in enumerate(constraints):
        for col, coeff in row.items():
            matrix[i, col] = coeff
    bounds = [(None, None)] + [(-1, 1)] * len(columns)
    ineq = np.c_[-np.ones(2 * len(matrix)), np.vstack([matrix, -matrix])]
    rhs = np.r_[-np.array(constants), constants]
    result = linprog(np.r_[1, np.zeros(len(columns))], A_ub=ineq, b_ub=rhs,
                     bounds=bounds, method="highs")
    assert result.success, result.message
    expected = sum(abs(b) * e * (p*lag/n + max(0, p-k)*(1-lag/n))
                   for n, lag, p, k, b, e in agents)
    assert abs(result.fun - expected) < 1e-9, (agents, result.fun, expected)
    return result.fun


def nonlinear_replay():
    # Integrate the original ODE across both information regimes. On each
    # constant-command interval the pitch polynomial is exact. Use Gaussian
    # quadrature for terminal displacement; H=A=B=1 in this test.
    nodes, weights = leggauss(12)
    value, K, _ = support(1, 1)
    n, lag = 4, 1
    for d in (0.02, K, 0.2):
        price = lag/n*K + (1-lag/n)*max(0, K-d)
        errors = []
        for e in (0.01, 0.001, 0.0001):
            h = 1-e
            q = one_switch(h, 1)
            worst = float("inf")
            for world in product((-d, 0.0, d), repeat=n):
                for sign in (-1, 1):
                    arcs = []
                    for j, wind in enumerate(world):
                        sensed = world[j-lag] if j >= lag else 0
                        diff = -np.clip(sensed/K, -1, 1)
                        arcs.append((e/n, 1-abs(diff)/2, diff, wind))
                    arcs.extend([(q, 0.5, sign, -sign*d), (h-q, 1, 0, -sign*d)])
                    t = pitch = rate = x_end = 0.0
                    for length, u, v, w in arcs:
                        local = (nodes+1)*length/2
                        angle = pitch + rate*local + v*local**2/2
                        x_end += length/2 * np.dot(weights, (1-t-local)*(u*np.sin(angle)+w))
                        pitch += rate*length + v*length**2/2
                        rate += v*length
                        t += length
                    worst = min(worst, sign*x_end)
            measured = (value-d/2-worst)/e
            errors.append(abs(measured-price))
        assert errors[-1] < 5e-5 and errors[-1] < errors[0]/20, errors
        print(f"nonlinear replay d={d:.8f} slope={measured:.8f} limit={price:.8f}")


def check_pulses():
    # Exact first two input moments for a discontinuous, non-binary reference.
    # Evaluate both trajectories on a denser mesh including non-knot times.
    rng = np.random.default_rng(8472)
    nodes, weights = leggauss(8)
    worst_ratio = 0.0
    for m in (2, 4, 8, 16):
        h, B = 1/m, 1.7
        subdivisions = 11
        edges = np.linspace(0, 1, subdivisions*m+1)
        banks = rng.uniform(size=(subdivisions*m, 2))
        pulses = []
        for start in range(0, len(banks), subdivisions):
            widths = np.diff(edges[start:start+subdivisions+1])
            centers = edges[start:start+subdivisions] + widths/2
            local = banks[start:start+subdivisions]
            mass = widths @ local
            center = (widths*centers) @ local / mass
            pulses.append((center-mass/2, center+mass/2))
        events = np.unique(np.r_[edges, np.asarray(pulses).ravel(),
                                 np.linspace(0, 1, 79)])
        angle = np.zeros(2)
        rate = np.zeros(2)
        pos = np.zeros(2)
        speed = np.zeros(2)
        for left, right in zip(events, events[1:]):
            mid = (left+right)/2
            relative = (nodes+1)*(right-left)/2
            pstart, pend = pulses[min(int(mid*m), m-1)]
            force = np.stack((banks[min(int(mid*m*subdivisions), len(banks)-1)],
                              ((pstart <= mid) & (mid < pend)).astype(float)))
            v = force[:, 1] - force[:, 0]
            thrust = force.mean(axis=1)
            theta = angle[:, None] + np.outer(rate, relative) + B*np.outer(v, relative**2)/2
            acceleration = thrust[:, None]*np.sin(theta)
            length = right-left
            pos += speed*length + length/2 * (acceleration @ (weights*(length-relative)))
            speed += length/2 * (acceleration @ weights)
            angle += rate*length + B*v*length**2/2
            rate += B*v*length
            if abs(right*m-round(right*m)) < 1e-12:
                assert abs(angle[1]-angle[0]) < 1e-12 and abs(rate[1]-rate[0]) < 1e-12
            bound = (1+4*B+B**2)*h**2/2
            worst_ratio = max(worst_ratio, abs(pos[1]-pos[0])/bound)
    assert 0 < worst_ratio < 1, worst_ratio
    print(f"pulse moment check passed; maximum error/bound={worst_ratio:.7f}")


def search_for_two_switch_counterexample():
    # n=2: compare the exact two checkpoints to every three-piece reduced
    # command with cuts on a 16-cell grid and amplitudes on a five-level grid.
    # |v| reduction preserves the switch restriction and improves both positions.
    split, B = 0.25, 1.0
    checkpoint, _, q = support(split, B)
    speed = quad(lambda t: (0.5 if t < q else 1.0) *
                 sin(B * (min(t, q)*t-min(t, q)**2/2)), 0, split, points=[q])[0]
    suffix, _, _ = support(1-split, B, B*(q*split-q**2/2), B*q)
    final = checkpoint + speed*(1-split) + suffix
    nodes, weights = leggauss(6)
    t = ((np.arange(16)[:, None] + (1+nodes)/2)/16).ravel()
    age = np.maximum(t[:, None]-np.arange(16)/16, 0)
    kernel = B*(age**2-np.maximum(age-1/16, 0)**2)/2
    weights = np.tile(weights/32, 16)
    objectives = np.c_[weights*np.maximum(split-t, 0), weights*(1-t)]
    amplitudes = np.array(list(product(np.linspace(0, 1, 5), repeat=3)))
    smallest = float("inf")
    for a, b in combinations(range(1, 16), 2):
        k = amplitudes[:, (np.arange(16) >= a).astype(int) + (np.arange(16) >= b)]
        acceleration = np.repeat(1-k/2, len(nodes), axis=1)*np.sin(k @ kernel.T)
        deficits = np.array([checkpoint, final]) - acceleration @ objectives
        smallest = min(smallest, np.min(np.max(deficits, axis=1)))
    assert smallest > 1e-9
    print(f"13125 two-switch competitors: minimum checkpoint deficit={smallest:.9g}")


if __name__ == "__main__":
    print("NumPy", np.__version__, "SciPy", scipy.__version__)
    check_rotors()
    for p in (0, 0.4, 1, 1.7):
        for lag in range(6):
            history_lp([(5, lag, p, 1, 1, 1)])
    print("24 single-aircraft minimax LPs passed")
    residual = history_lp([(4, 1, 1.5, 1, -1, 1), (3, 2, 0.4, 1, 1, 0.6)])
    print(f"two-aircraft, 128 joint histories, shared-sign LP residual={residual:.9f}")
    nonlinear_replay()
    check_pulses()
    search_for_two_switch_counterexample()
