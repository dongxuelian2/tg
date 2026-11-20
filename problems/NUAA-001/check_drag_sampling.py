"""Compare opposite, zero-report winds to the continuous L1 estimation bound."""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import root


def check(h, polynomial=False):
    # Normalize the first two exponential kernels into an exact stable basis
    # for rates 1,2: 1-exp(-2s)=2(1-exp(-s))-(1-exp(-s))**2.
    def basis(x):
        y = x if polynomial else -np.expm1(-h*x)/h
        return np.array([y, y*y])

    def integrate(f, x=0., y=1.):
        return quad(f, x, y, epsabs=1e-13)[0]

    def cancel(cuts):
        a, b = cuts
        return [integrate(lambda x: basis(x)[k]) -
                2*integrate(lambda x: basis(x)[k], a, b) for k in range(2)]

    sol = root(cancel, [.45, .85], tol=1e-12)
    a, b = sol.x
    assert 0 < a < b < 1 and max(abs(np.array(cancel(sol.x)))) < 1e-11
    coeff = np.linalg.solve(np.stack([basis(a), basis(b)]), np.ones(2))
    primal = sum(integrate(lambda x: abs(1-coeff @ basis(x)), lo, hi)
                 for lo, hi in [(0,a), (a,b), (b,1)])
    dual = 1-2*(b-a)
    assert abs(primal-dual) < 1e-10
    return dual


if __name__ == "__main__":
    limit = check(1., polynomial=True)
    assert 1/9 < limit < 1/4
    for h in [1., .1, .01, .001]:
        print(h, check(h))
    print("two-rate normalized short-window limit:", limit)
