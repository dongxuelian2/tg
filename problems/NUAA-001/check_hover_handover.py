"""Compare the Chebyshev handover bound to moment LPs and unscaled rigid-body ODEs."""

import numpy as np
from scipy.optimize import linprog
from check_full_body import decay, e3, endpoint, hover, inertia, initial, mixer


def moments(nodes):
    return np.array([np.diff(nodes**k)/k for k in range(1, 6)])


if __name__ == "__main__":
    switches = (1-np.cos(np.arange(6)*np.pi/5))/2
    signs = np.array([1, -1, 1, -1, 1])
    vals = moments(switches) @ signs
    assert np.max(np.abs(vals[:4])) < 1e-14 and abs(vals[4]-1/256) < 1e-14
    for count in (32, 128, 512):
        mat = moments(np.linspace(0, 1, count+1))
        lp = linprog(-mat[4], A_eq=mat[:4], b_eq=np.zeros(4), bounds=(-1, 1))
        assert lp.success and abs(-lp.fun-1/256) < .07/count**2
    q = np.array([-.12, -.12, .12, .12])
    target = np.cross((mixer[1:]@(decay*q))/inertia, e3)[:2]/(24*256)
    errors = []
    for h in (.16, .08, .04):
        state = initial.copy()
        for ds, sign in zip(np.diff(switches), signs):
            state = endpoint(h*ds, lambda s: hover+sign*q, state)
        errors.append(np.linalg.norm(state[:2]/h**5-target))
    assert all(b < .6*a for a, b in zip(errors, errors[1:])), errors
    print("Endpoint moment checks passed. Uncorrected scaled displacement errors:", errors)
