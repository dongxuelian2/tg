"""Numerical counterexample searches for the finite-motor-delay working lemma."""

from itertools import product
from math import cos, exp, expm1, sin

import numpy as np
from scipy.integrate import quad
from scipy.optimize import linprog, minimize_scalar


def support(h, D, hover, delta_sum=0.0, delta_diff=0.0):
    def angle_kernel(t):
        return t*t/2-t/D-expm1(-D*t)/D**2

    def evaluate(q, coefficients=False):
        def integrand(t):
            residual = exp(-D*t)
            u = 0.5+(hover+delta_sum-0.5)*residual
            u += -0.5*expm1(-D*max(t-q, 0))
            impulse = t/D+expm1(-D*t)/D**2
            angle = angle_kernel(t)-angle_kernel(max(t-q, 0))+delta_diff*impulse
            if coefficients:
                return ((h-t)*D*residual*sin(angle), (h-t)*D*u*cos(angle)*impulse)
            return (h-t)*u*sin(angle)
        if coefficients:
            return [quad(lambda t: integrand(t)[i], 0, h, points=[q])[0] for i in (0, 1)]
        return quad(integrand, 0, h, points=[q])[0]

    best = minimize_scalar(lambda q: -evaluate(q), bounds=(0, h), method="bounded",
                           options={"xatol": 3e-14})
    assert best.success and 0 < best.x < h
    return -best.fun, *evaluate(best.x, True), best.x


def joint_linear_program(agents):
    # P, k0, gamma, slots, feedback delay in slots, signed weight
    histories = [list(product((-1, 1), repeat=n)) for _, _, _, n, _, _ in agents]
    variables = {}
    rows = []
    for joint in product(*histories):
        wind = 0
        choices = []
        for i, ((P, k0, gamma, n, lag, b), word) in enumerate(zip(agents, joint)):
            wind += b*P*sum(word)/n
            for t in range(n):
                key = (i, t, word[:max(0, t+1-lag)])
                j = variables.setdefault(key, len(variables))
                choices.append((j, b*k0/n, abs(b)*gamma/n))
        for sign in (-1, 1):
            rows.append((choices, sign, sign*wind))
    count = len(variables)
    matrix = np.zeros((len(rows)+2*count, 2*count+1))
    rhs = np.zeros(len(matrix))
    for r, (choices, sign, bound) in enumerate(rows):
        matrix[r, -1] = -1
        rhs[r] = bound
        for j, gain, cost in choices:
            matrix[r, j] = -sign*gain
            matrix[r, count+j] = cost
    for j in range(count):
        r = len(rows)+2*j
        matrix[r:r+2, j] = [1, -1]
        matrix[r:r+2, count+j] = -1
    result = linprog(np.r_[np.zeros(2*count), 1], A_ub=matrix, b_ub=rhs,
                     bounds=[(-1, 1)]*count+[(0, 1)]*count+[(None, None)])
    expected = sum(abs(b)*(P*lag/n+(1-lag/n)*(max(0, P-k0)+gamma*min(P/k0, 1)))
                   for P, k0, gamma, n, lag, b in agents)
    assert result.success and abs(result.fun-expected) < 1e-8, (result, expected)


if __name__ == "__main__":
    for D in (0.4, 2.0, 20.0):
        for hover in (0.25, 0.7):
            base, c, k0, q = support(1, D, hover)
            eps = 1e-5
            numeric = (support(1+eps, D, hover)[0]-support(1-eps, D, hover)[0])/(2*eps)
            assert abs(numeric-(k0+c*(0.5-hover))) < 1e-7 and k0 > c/2 > 0
            for i, expected in enumerate((c, k0)):
                vector = np.eye(2)[i]*eps
                central = (support(1, D, hover, *vector)[0]-
                           support(1, D, hover, *(-vector))[0])/(2*eps)
                assert abs(D*central-expected) < 1e-7
            print(f"D={D:g} u0={hover:g} q={q:.7f} k0={k0:.7f} gamma={c/2:.7f}")
    for P, lag in product((0, 0.3, 1, 2.1), range(5)):
        joint_linear_program([(P, 1, 0.22, 4, lag, 1)])
    joint_linear_program([(0.35, 1, 0.21, 3, 1, -1.4),
                          (2.0, 1.1, 0.6, 4, 3, 0.8)])
    print("20 single-aircraft LPs and one 128-history joint minimax LP passed")
