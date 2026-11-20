"""Finite-time 16-component shooting in the unscaled quaternion/motor model."""

import numpy as np
from scipy.linalg import expm
from scipy.optimize import linprog
from check_full_body import (decay, e3, endpoint, hover, hover_linearization,
                             inertia, initial, mixer, scaled)


def verify(h):
    breaks = (1-np.cos(np.arange(6)*np.pi/5))/2
    durations = h*np.diff(breaks)
    # Include collective and yaw components and off-axis horizontal translation.
    q = np.array([-.095, .14, .085, -.04])
    commands = hover + np.array([1,-1,1,-1,1])[:, None]*q
    target = np.zeros(16)
    target[:2] = np.cross((mixer[1:]@(decay*q))/inertia, e3)[:2]/(24*256)

    def forward(table):
        state = initial.copy()
        for dt, cmd in zip(durations, table):
            state = endpoint(dt, lambda s: cmd, state)
        return scaled(state, h)

    # Linearized finite-time impulse map at hover, not fitted to the solution.
    A, B = hover_linearization(1)
    aug = np.zeros((20, 20))
    aug[:16, :16], aug[:16, 16:] = A, B
    blocks = [expm(A*h*(1-breaks[i+1])) @ expm(aug*dt)[:16, 16:]
              for i, dt in enumerate(durations)]
    powers = [5,5,3,4,4,2,3,3,3,2,2,2,1,1,1,1]
    jac = np.hstack(blocks)/h**np.array(powers)[:, None]
    start = np.max(np.abs(forward(commands)-target))
    for step in range(10):
        residual = target-forward(commands)
        if max(abs(residual)) < 2e-10:
            break
        # Minimize the largest correction, while enforcing the actual command box.
        bounds = [(-c, .5-c) for c in commands.ravel()]+[(0, None)]
        solution = linprog([0]*20+[1], A_eq=np.c_[jac, np.zeros(16)], b_eq=residual,
                           A_ub=np.c_[np.r_[np.eye(20), -np.eye(20)], -np.ones(40)],
                           b_ub=np.zeros(40), bounds=bounds,
                           options={"primal_feasibility_tolerance":1e-10})
        assert solution.success, solution.message
        commands += solution.x[:20].reshape(5, 4)
    remaining = max(abs(target-forward(commands)))
    margin = min(commands.min(), .5-commands.max())
    assert remaining < 2e-10 and margin > .05
    return h, start, remaining, margin


if __name__ == "__main__":
    for duration in (.12, .24, .48):
        print("H, initial/final scaled residual, command headroom:", verify(duration))
