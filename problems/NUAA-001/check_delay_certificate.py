"""Rational interval verification of a three-pulse delayed-state safety policy."""
from fractions import Fraction as Q
from math import factorial


def add(x,y):
    return x[0]+y[0], x[1]+y[1]


def mul(x,y):
    ends=[a*b for a in x for b in y]
    return min(ends),max(ends)


def point(x):
    return Q(x),Q(x)


def neg(x):
    return -x[1],-x[0]


def mag(x):
    return max(abs(x[0]),abs(x[1]))


def verify():
    one,zero=point(1),point(0)
    upper=sum(Q((-1)**k,factorial(k)) for k in range(17))
    decay=(upper-Q(1,factorial(17)),upper)
    complement=add(one,neg(decay))
    F=[[one,one,decay],[zero,one,complement],[zero,zero,decay]]
    B=[add(point(Q(1,2)),neg(decay)),decay,complement]
    gain=[[-1582,-3034,-1649],[2164,2569,799],[-582,-534,-149]]
    P=[Q(7,1000),Q(4,1000),Q(1,1000)]
    M=[[point(int(i==j)) for j in range(3)] for i in range(3)]
    max_command=0
    for row in gain:
        u=[Q(c,1000) for c in row]
        max_command=max(max_command,sum(abs(c)*r for c,r in zip(u,P)))
        M=[[add(sum_intervals(mul(F[i][k],M[k][j]) for k in range(3)),
                mul(B[i],point(u[j]))) for j in range(3)] for i in range(3)]
    nominal=[sum(mag(m)*r for m,r in zip(row,P)) for row in M]
    worst_wind=[Q(6,1000),Q(3,1000),Q(0)]
    contraction=max(x/r for x,r in zip(nominal,P))
    assert all(x+y<r for x,y,r in zip(nominal,worst_wind,P))
    assert contraction<Q(1,100) and max_command<Q(27,1000)
    # The actual filtered acceleration is bounded by max(max|a0|, max|u|).
    intersample=P[0]+3*P[1]+Q(9,2)*max(P[2],max_command)+Q(1,2000)*Q(7,2)**2
    assert intersample<Q(16,100)
    print("certified upper bounds: command, intersample height, contraction")
    print(*(float(x) for x in [max_command,intersample,contraction]))


def sum_intervals(terms):
    result=point(0)
    for term in terms:
        result=add(result,term)
    return result


if __name__ == "__main__":
    verify()
