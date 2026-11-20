"""Falsification checks for the local NUAA-001 infinite-horizon working notes."""

from itertools import product
from math import exp

import numpy as np
from scipy.spatial.transform import Rotation


def skew(x):
    a, b, c = x
    return np.array([[0, -c, b], [c, 0, -a], [-b, a, 0]])


def right_jacobian(a):
    radius = np.linalg.norm(a)
    cross = skew(a)
    if radius < 1e-4:
        c1, c2 = 0.5-radius**2/24, 1/6-radius**2/120
    else:
        c1 = (1-np.cos(radius))/radius**2
        c2 = (radius-np.sin(radius))/radius**3
    return np.eye(3)-c1*cross+c2*cross@cross


def hermite(h, t, a, b, u, v):
    s = t/h
    position = ((1+2*s)*(1-s)**2*a + s**2*(3-2*s)*b
                + h*s*(1-s)**2*u + h*s**2*(s-1)*v)
    speed = ((6*s-6*s*s)*(b-a)/h + (1-4*s+3*s*s)*u + (3*s*s-2*s)*v)
    acceleration = ((6-12*s)*(b-a)/h**2 + (6*s-4)*u/h + (6*s-2)*v/h)
    return position, speed, acceleration


def geometry_tests():
    rng = np.random.default_rng(1729)
    largest_residual = 0
    for h in (0.02, 0.3, 1.0):
        for _ in range(8):
            omega = rng.normal(size=3)
            endpoints = rng.uniform(-0.02, 0.02, (2, 3))
            # Endpoints specify the new body rates in their own body coordinates.
            rates = omega+rng.uniform(-0.03, 0.03, (2, 3))
            b = [np.linalg.solve(right_jacobian(a), om-Rotation.from_rotvec(a).inv().apply(omega))
                 for a, om in zip(endpoints, rates)]
            base = Rotation.random(random_state=rng)

            def output(t):
                a, da, dda = hermite(h, t, *endpoints, *b)
                adjustment = Rotation.from_rotvec(a)
                pose = base * Rotation.from_rotvec(t*omega) * adjustment
                rate = adjustment.inv().apply(omega) + right_jacobian(a)@da
                return pose, rate, a, da, dda

            for j in range(2):
                pose, rate, *_ = output(j*h)
                expected = base*Rotation.from_rotvec(j*h*omega)*Rotation.from_rotvec(endpoints[j])
                assert (expected.inv()*pose).magnitude() < 1e-12
                assert np.linalg.norm(rate-rates[j]) < 1e-11
            A = max(np.linalg.norm(a) for a in endpoints)
            B = max(np.linalg.norm(db) for db in b)
            step = 1e-5*h
            for t in np.linspace(h/100, 99*h/100, 29):
                _, rate, a, da, dda = output(t)
                relative = output(t-step)[0].inv()*output(t+step)[0]
                numerical = relative.as_rotvec()/(2*step)
                largest_residual = max(largest_residual, np.linalg.norm(numerical-rate))
                assert np.linalg.norm(a) <= A+2*h*B+1e-11
                assert np.linalg.norm(da) <= 3*A/h+2*B+1e-11
                assert np.linalg.norm(dda) <= 12*A/h**2+8*B/h+1e-11
    assert largest_residual < 1e-7
    print(f"SO(3) endpoint and intersample checks passed; max rate residual {largest_residual:.3g}")


def motor_test():
    # Scalar extremal recurrence saturates the motor bound.
    for D, h, jump, command in product((0.2, 7.), (.01, .4), (.0001, .02), (0., .01)):
        rho = exp(-D*h)
        start = .003
        err = start
        for j in range(1, 401):
            err = rho*err + (1-rho)*command + jump
            bound = rho**j*start+(1-rho**j)*command+jump*(1-rho**j)/(1-rho)
            assert abs(err-bound) < 1e-11
    print("motor reset bounds checked including worst-sign accumulation")


def relation(a, b):
    return (b[0]-a[0], b[1]-a[1])


def communication_tests():
    moves = ((0, 0), (1, 0), (0, 1))
    count = 0
    for n in range(1, 7):
        for a in product(range(2**n), repeat=2):
            for j, delta in enumerate(moves):
                b = tuple(ax+dx for ax, dx in zip(a, delta))
                if max(b) >= 2**n:
                    continue
                # Generate the stream from input digits and a carry, without
                # consulting the two final coordinates for the relation.
                carry = list(delta)
                residues = [0, 0]
                for k in range(n):
                    digits = [(ax >> k) & 1 for ax in a]
                    r = [digits[i]+carry[i] for i in range(2)]
                    carry = [x//2 for x in r]
                    out = [x % 2 for x in r]
                    assert out == [(x >> k) & 1 for x in b]
                    for i, pair in enumerate((digits, out)):
                        residues[i] += pow(2, k, 3)*(pair[0]+2*pair[1])
                assert not any(carry) and (residues[1]-residues[0]) % 3 == j
                assert relation(a, b) == moves[j]
                count += 1
    print(f"private-coordinate monitor and ternary disclosure verified on {count} word pairs")


if __name__ == "__main__":
    geometry_tests()
    motor_test()
    communication_tests()
