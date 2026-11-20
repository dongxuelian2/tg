"""Finite cross-checks for permanent-synthesis.md; run directly with Python 3."""

from itertools import product, zip_longest
from math import exp, expm1

import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation
import z3


def coordinate_word(x, t):
    return [(int(j <= x), int(j <= t)) for j in range(1, max(x, t) + 1)]


def monitor(left, right):
    # One finite flag per possible relation; no coordinate counters.
    candidates = {(0, 0), (1, 0), (0, 1)}
    extra = [0, 0]
    for a, b in zip_longest(left, right, fillvalue=(0, 0)):
        for k in (0, 1):
            if b[k] < a[k]:
                return None
            extra[k] = min(2, extra[k] + b[k] - a[k])
        candidates = {r for r in candidates if all(extra[k] <= r[k] for k in (0, 1))}
    return tuple(extra) if tuple(extra) in candidates else None


def check_queries():
    coordinates = list(product(range(9), repeat=2))
    trials = 0
    for left, right in product(coordinates, repeat=2):
        a, b = coordinate_word(*left), coordinate_word(*right)
        diff = tuple(y - x for x, y in zip(left, right))
        expected = diff if diff in ((0, 0), (1, 0), (0, 1)) else None
        assert monitor(a, b) == expected
        if expected is not None:
            end_a, end_b = 3 * (len(a) + 1), 3 * (len(b) + 1)
            assert max(end_a, end_b) < min(end_a, end_b) + 4
        trials += 1
    return trials


def check_tableau(width=5):
    # A two-state head scans right and rejects if it reads consecutive ones.
    # Encode a cell as 2*head + bit, with head 0 (absent), 1, 2 or 3 (reject).
    # The solver sees overlapping triples, not the high-level string condition.
    rows = [[[z3.Int(f"tile_{x}_{t}_{i}") for i in range(3)]
             for x in range(width)] for t in range(width + 1)]
    solver = z3.Solver()
    for t, row in enumerate(rows):
        for x, (left, center, right) in enumerate(row):
            solver.add(center >= 0, center < 6)
            if x == 0:
                solver.add(left == -2)
            if x == width - 1:
                solver.add(right == 0)
            else:
                solver.add(center == row[x + 1][0], right == row[x + 1][1])
            if t == width:
                continue
            reject = z3.And(center / 2 == 2, center % 2 == 1)
            arrives = z3.Or(left / 2 == 1, z3.And(left / 2 == 2, left % 2 == 0))
            new_head = z3.If(reject, 3, z3.If(arrives, 1 + left % 2, 0))
            solver.add(rows[t + 1][x][1] == 2 * new_head + center % 2)
    count = 0
    for bits in product(range(2), repeat=width):
        solver.push()
        for x, bit in enumerate(bits):
            solver.add(rows[0][x][1] == (2 if x == 0 else 0) + bit)
        status = solver.check()
        assert status in (z3.sat, z3.unsat)
        expected = all(not (a and b) for a, b in zip(bits, bits[1:]))
        assert (status == z3.sat) == expected
        count += int(expected)
        solver.pop()
    return count


def encoders(size):
    # Each table indexes (private alarm, private data), in that order.
    for code in product(range(size), repeat=4):
        if set(code[:2]).isdisjoint(code[2:]):
            yield code


def parity_satisfiable(equations):
    # Variables are local decoder entries. Propagate prescribed XOR on edges.
    neighbors = {}
    for u, v, parity in equations:
        neighbors.setdefault(u, []).append((v, parity))
        neighbors.setdefault(v, []).append((u, parity))
    values = {}
    for root in neighbors:
        if root in values:
            continue
        values[root] = 0
        stack = [root]
        while stack:
            node = stack.pop()
            for other, parity in neighbors[node]:
                required = values[node] ^ parity
                if other not in values:
                    values[other] = required
                    stack.append(other)
                elif values[other] != required:
                    return False
    return True


def has_protocol(size, target):
    codes = list(encoders(size))
    for emit_a, emit_b in product(codes, repeat=2):
        equations = []
        for a, b, x, y in product(range(2), repeat=4):
            # Local sent letter can be reconstructed from the same local input.
            u = (0, a, x, emit_b[2 * b + y])
            v = (1, b, y, emit_a[2 * a + x])
            equations.append((u, v, target(x, y)))
        if parity_satisfiable(equations):
            return True
    return False


def check_capacity():
    for k in (1, 2, 3, 4):
        assert has_protocol(k, lambda x, y: x & y) == (k >= 4)
        assert has_protocol(k, lambda x, y: x ^ y) == (k >= 2)


def hermite_acceleration(s, h, start, velocity):
    # Relative position/velocity start at these values and both vanish at h.
    r = s / h
    return (12 * r - 6) * start / h**2 + (6 * r - 4) * velocity / h


def hermite_value(s, h, start, velocity):
    r = s / h
    return (1 - r)**2 * ((1 + 2 * r) * start + s * velocity)


def thrust_direction(angle):
    return np.array([np.sin(angle), np.cos(angle)])


def check_motor_shadow(eta):
    # A nonlinear pitch manoeuvre with four independent first-order motor states.
    # State order is horizontal/vertical position, velocity, pitch, rate, thrusts.
    h, rate = 0.4, 3.0
    real = np.r_[0.0032, -0.0015, 0.013, -0.014, 0.023, 0.042, np.ones(4)]
    nominal = real.copy()
    max_endpoint_error = max_wrench = 0.0
    for frame in range(12):
        nominal = eta * np.round(nominal / eta)
        initial_gap = real - nominal
        sign = 1 if frame % 2 else -1
        command = np.array([1 + 0.22 * sign, 1.03, 1 - 0.19 * sign, 1.03])

        def system(s, both):
            nonlocal max_wrench
            n, a = np.split(both, 2)
            # Compute the imposed wrench from nominal data and start error only.
            angle_offset = hermite_value(s, h, initial_gap[4], initial_gap[5])
            angle_accel = hermite_acceleration(s, h, initial_gap[4], initial_gap[5])
            motor_offset = initial_gap[6:] * exp(-rate * s)
            accel_offset = hermite_acceleration(s, h, initial_gap[:2], initial_gap[2:4])
            thrust_n = n[6:].sum()
            force = (accel_offset + thrust_n * thrust_direction(n[4])
                     - (thrust_n + motor_offset.sum()) * thrust_direction(n[4] + angle_offset))
            torque = angle_accel - motor_offset[0] + motor_offset[2]
            max_wrench = max(max_wrench, np.linalg.norm(force), abs(torque))

            def derivative(state, external_force, external_torque):
                thrusts = state[6:]
                acceleration = thrusts.sum() * thrust_direction(state[4]) - [0, 4]
                return np.r_[state[2:4], acceleration + external_force, state[5],
                             thrusts[0] - thrusts[2] + external_torque,
                             rate * (command - thrusts)]

            return np.r_[derivative(n, 0, 0), derivative(a, force, torque)]

        solution = solve_ivp(system, (0, h), np.r_[nominal, real], method="DOP853",
                             atol=1e-13, rtol=1e-12)
        assert solution.success
        nominal, real = np.split(solution.y[:, -1], 2)
        final_gap = real - nominal
        max_endpoint_error = max(max_endpoint_error, np.max(np.abs(final_gap[:6])))
        np.testing.assert_allclose(final_gap[6:], initial_gap[6:] * exp(-rate * h),
                                   atol=1e-12, rtol=0)
        assert np.max(np.abs(initial_gap[6:])) <= eta / (-expm1(-rate * h)) + 1e-12
    assert max_endpoint_error < 1e-10 and max_wrench < 100 * eta
    return max_endpoint_error, max_wrench


def check_shadow():
    measurements = [(step, *check_motor_shadow(step)) for step in (0.003, 0.0003)]
    assert measurements[1][2] < measurements[0][2] / 4
    return measurements


def check_attitude_endpoints():
    h = 0.7
    spin = np.array([0.31, -0.24, 0.43])
    start = np.array([0.08, -0.05, 0.03])
    actual_spin = spin + [0.011, 0.024, -0.019]
    skew = np.cross(start, np.eye(3)).T
    angle = np.linalg.norm(start)
    jacobian = (np.eye(3) - (1 - np.cos(angle)) / angle**2 * skew
                + (angle - np.sin(angle)) / angle**3 * (skew @ skew))
    initial_rate = np.linalg.solve(jacobian,
                                   actual_spin - Rotation.from_rotvec(-start).apply(spin))

    def attitude(t):
        return (Rotation.from_rotvec(t * spin)
                * Rotation.from_rotvec(hermite_value(t, h, start, initial_rate)))

    for t, expected in ((0, actual_spin), (h, spin)):
        step = 1e-6
        increment = attitude(t - step).inv() * attitude(t + step)
        np.testing.assert_allclose(increment.as_rotvec() / (2 * step), expected,
                                   atol=2e-9, rtol=0)


if __name__ == "__main__":
    print("coordinate word pairs:", check_queries())
    print("accepted initial rows of the five-cell tableau:", check_tableau())
    check_capacity()
    print("one-frame encoder and decoder checks passed")
    check_attitude_endpoints()
    print("noncommuting attitude endpoint rates passed")
    for step, error, wrench in check_shadow():
        print(f"grid={step:g}, mechanical endpoint error={error:.3g}, "
              f"largest sampled force/torque norm={wrench:.4g}")
