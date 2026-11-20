"""Independently integrate the actuator and compare window certificates to LP."""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linprog


def travel(h, rate):
    return h/rate * (np.logaddexp(rate*h/4, -rate*h/4)-np.log(2))


def check_transfer(h, rate):
    switch = (np.logaddexp(0, rate*h/2)-np.log(2))/rate
    state = np.zeros(3)
    for t0, t1, command in [(0, switch, 1), (switch, switch+h/2, -1),
                            (switch+h/2, h, 1)]:
        def rhs(t, x):
            return [x[1], x[2], rate*(command-x[2])]
        result = solve_ivp(rhs, (t0, t1), state, atol=1e-13, rtol=1e-11)
        assert result.success
        state = result.y[:, -1]
    assert np.max(np.abs(state-[travel(h, rate), 0, 0])) < 1e-9


def check_windows(rng):
    limits = travel(rng.uniform(.5, 5, 4), 1.4)
    lower = rng.uniform(-2, 0.5, (2, 5))
    upper = lower + rng.uniform(0.5, 3.5, (2, 5))
    position = np.r_[0, np.cumsum(limits)]
    distance = np.abs(position[:, None] - position)
    heights = np.max(lower.max(axis=0)[:, None]-distance, axis=0)
    expected = np.all(heights <= upper.min(axis=0))
    adjacent = np.diff(np.eye(5), axis=0)
    result = linprog(np.zeros(5), A_ub=np.vstack([adjacent, -adjacent]),
                     b_ub=np.tile(limits, 2),
                     bounds=list(zip(lower.max(axis=0), upper.min(axis=0))))
    assert result.success == expected
    if expected:
        assert np.all(np.abs(np.diff(heights)) <= limits+1e-12)
    return bool(expected)


if __name__ == "__main__":
    for h in (.2, 1., 3., 7.):
        for rate in (.3, 1., 4.):
            check_transfer(h, rate)
    rng = np.random.default_rng(711)
    wins = sum(check_windows(rng) for _ in range(250))
    print(f"12 actuator integrations and 250 LPs passed ({wins} feasible)")
