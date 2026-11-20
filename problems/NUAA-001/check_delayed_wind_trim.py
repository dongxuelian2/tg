"""Grid primal and continuous-polynomial dual for the wind-trim moment problem."""

import numpy as np
from numpy.polynomial import Polynomial
from scipy.optimize import linprog


def bracket(fraction, cells=512):
    edges = np.linspace(0, 1, cells+1)
    moments = np.array([np.diff(edges**(j+1))/(j+1) for j in range(5)])
    coeff = np.eye(5)[:4].copy()
    coeff[2, 4] = -fraction**2/6
    coeff[3, 4] = -fraction/2
    constraints = coeff @ moments
    result = linprog(-moments[4], A_eq=constraints, b_eq=np.zeros(4),
                     bounds=(-1, 1), options={"dual_feasibility_tolerance": 1e-10})
    assert result.success and max(abs(constraints @ result.x)) < 1e-9
    residual = Polynomial(np.eye(5)[4] + result.eqlin.marginals @ coeff)
    roots = [r.real for r in residual.roots() if abs(r.imag)<1e-9 and 0<r.real<1]
    integral = residual.integ()
    cuts = [0, *sorted(roots), 1]
    upper = sum(abs(integral(b)-integral(a)) for a, b in zip(cuts, cuts[1:]))
    lower = -result.fun
    assert lower-1e-10 <= upper <= lower*1.001
    return lower/12, upper/12


if __name__ == "__main__":
    for fraction in (1., .5, .1, .01, .001):
        print("remaining fraction, numerical bracket for beta/L^5:",
              fraction, bracket(fraction))
