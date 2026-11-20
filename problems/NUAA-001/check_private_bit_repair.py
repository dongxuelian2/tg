"""Exact subset tests, integer-flow witness, and direct private-control LPs."""
from fractions import Fraction
from math import expm1

import numpy as np
from scipy.optimize import brentq, linprog
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow


edges = [(i, 5+(i+j) % 5) for i in range(5) for j in range(3)]
goals = np.r_[np.ones(5), -np.ones(5)]


def split(i, j, mask):
    return ((mask >> i) ^ (mask >> (j-5))) & 1


def threshold(mask):
    need = Fraction(0)
    for group in range(1, 1 << 10):
        boundary = inside = 0
        for i, j in edges:
            u, v = (group >> i) & 1, (group >> j) & 1
            boundary += u ^ v
            inside += u*v*split(i, j, mask)
        deficit = group.bit_count()-boundary
        if deficit > 0:
            if inside == 0:
                return None
            need = max(need, Fraction(deficit, 2*inside))
    return need


def direct(mask, late):
    cols, caps = [], []
    for i, j in edges:
        col = np.zeros(10)
        col[[i,j]] = 1
        cols.append(col)
        if split(i, j, mask):
            caps.append(1-late)
            for endpoint in (i, j):
                single = np.zeros(10)
                single[endpoint] = 1
                cols.append(single)
                caps.append(late)
        else:
            caps.append(1.)
    return linprog(np.zeros(len(cols)), A_eq=np.array(cols).T, b_eq=goals,
                   bounds=[(-k,k) for k in caps]).success


def integer_witness(mask):
    capacity = np.zeros((12, 12), dtype=int)
    capacity[10, :10] = 6
    for i, j in edges:
        local = 5*split(i, j, mask)
        capacity[i, j] = capacity[j, i] = 6-local
        capacity[[i,j], 11] += local
    flow = maximum_flow(csr_matrix(capacity), 10, 11)
    witness = flow.flow.toarray()
    assert flow.flow_value == 60 and np.all(witness <= capacity)
    assert np.all(witness[:10].sum(axis=1) == 0)
    return flow.flow_value


if __name__ == "__main__":
    finite = []
    for mask in range(32):
        exact = threshold(mask)
        for late in (0., .5, 5/6, .95, 1.):
            assert direct(mask, late) == (exact is not None and float(exact) <= late)
        if exact is not None:
            finite.append(exact)
    assert min(finite) == threshold(3) == Fraction(5, 6)
    def impulse(s):
        return s*s/2-s-expm1(-s)
    deadline = brentq(lambda t: impulse(1-t)/impulse(1)-5/6, 0, 1)
    print("32 codes, 160 LP checks; optimum 5/6; flow", integer_witness(3),
          "of 60; relative deadline", deadline)
