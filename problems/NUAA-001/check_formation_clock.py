"""Integrate the planar invariant subsystem and check the motor-clock formulas."""

import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate import solve_ivp


clock = Polynomial([0.4, 0.8, 0.2, -0.12, 0.03, 0.02, -0.006])
jets = [clock.deriv(k) for k in range(6)]


def reference(t, radius):
    x, velocity, a, j, k, ell = [radius * jet(t) for jet in jets]
    d = 1 + a*a
    pitch = np.arctan(a)
    rate = j/d
    torque = k/d - 2*a*j*j/d**2
    torque_rate = ell/d - 6*a*j*k/d**2 + (6*a*a-2)*j**3/d**3
    thrust = np.sqrt(d)
    pair_sign = np.array([1., -1.])
    motors = (thrust + pair_sign*torque)/4
    commands = motors + (a*j/thrust + pair_sign*torque_rate)/4
    return np.array([x, velocity, 0, 0, pitch, rate, *motors]), commands


def simulate(radius):
    def dynamics(t, state):
        _, vx, _, vz, pitch, rate, plus, minus = state
        thrust = 2*(plus + minus)
        _, command = reference(t, radius)
        return [vx, thrust*np.sin(pitch), vz, thrust*np.cos(pitch)-1,
                rate, 2*(plus-minus), command[0]-plus, command[1]-minus]

    start, finish = -0.15, 0.25
    times = np.linspace(start, finish, 81)
    result = solve_ivp(dynamics, (start, finish), reference(start, radius)[0],
                       t_eval=times, rtol=1e-11, atol=1e-13, method="DOP853")
    assert result.success
    expected = np.stack([reference(t, radius)[0] for t in times]).T
    return np.max(np.abs(result.y - expected))


if __name__ == "__main__":
    radii = np.array([0., 1/np.sqrt(2), 2.])
    errors = [simulate(r) for r in radii]
    assert max(errors) < 1e-9
    thrust = np.sqrt(1 + radii**2)
    amplitude = 3.5*radii/(1+radii**2)
    low, high = (thrust-amplitude)/4, (thrust+amplitude)/4
    assert low[1] < 0 < min(low[0], low[2]) and max(high) < 1
    fine = np.linspace(0.01, 2, 100001)
    bounds = (1+fine**2)**1.5/fine
    assert abs(bounds.min() - 3*np.sqrt(3)/2) < 1e-8
    print("Maximum state integration error:", max(errors))
    for r, lower, upper in zip(radii, low, high):
        print(f"r={r:.8f}  commands=[{lower:.8f}, {upper:.8f}]")
