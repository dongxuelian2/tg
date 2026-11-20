"""Check stopping envelopes against direct ODE integration and a command LP."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, linprog


def peak(v, a, b=0.0):
    # Dimensionless D=U=1. Velocity under full braking has at most one maximum.
    k = 1-b
    vel = lambda t: v-k*t-(a+1)*np.expm1(-t)
    turn = max(0.0, np.log((a+1)/k)) if a > -1 else 0.0
    if vel(turn) <= 0:
        return 0.0, 0.0
    end = turn+1
    while vel(end) > 0:
        end *= 2
    t = brentq(vel, turn, end)
    displacement = v*t-k*t*t/2+(a+1)*(t+np.expm1(-t))
    return (displacement, t) if displacement > 0 else (0.0, 0.0)


def verify():
    tstar = brentq(lambda t: 2*(-np.expm1(-t))-t, 1, 2)
    width = tstar-tstar*tstar/2
    assert abs(peak(0., 1.)[0]-width) < 1e-12
    print("critical width / (U/D^2):", width)
    rng = np.random.default_rng(2891)
    worst = 0.
    for b in [0., .17, .55]:
        for _ in range(24):
            v,a = rng.uniform(-.8,.8),rng.uniform(-1,1)
            s,t = peak(v,a,b)
            # Separate event-based integration of the original motor and mechanical ODE.
            def at_peak(time, state):
                return state[1]
            at_peak.direction = -1
            ode = solve_ivp(lambda time,y: [y[1], y[2]+b, -1-y[2]],(0,12),
                            [0,v,a],events=at_peak,rtol=1e-11,atol=1e-13,max_step=.05)
            xs = [row[0] for row in ode.y_events[0]]
            worst=max(worst,abs(max([0]+xs)-s))
    assert worst < 1e-9
    print("72 ODE envelope checks, maximum error:",worst)

    # A witness for H=0.20 where each one-wall test is tight.
    H,v=.2,-.025
    a=brentq(lambda a: peak(v,a)[0]+peak(-v,-a)[0]-H,0,1)
    x,tm=peak(-v,-a)
    _,tp=peak(v,a)
    # For bounded controls independently optimized at every time, maximize the
    # minimum of the two wall clearances. The LP below tests a finite subfamily.
    knots=np.unique(np.r_[np.linspace(0,tp,101),tm])
    def row(t):
        kernel=lambda r:r*r/2-r-np.expm1(-r)
        rr=np.maximum(0,t-knots)
        return kernel(rr[:-1])-kernel(rr[1:])
    free=lambda t:x+v*t+a*(t+np.expm1(-t))
    matrix=np.array([np.r_[row(tp),1],np.r_[-row(tm),1]])
    res=linprog(np.r_[np.zeros(len(knots)-1),-1],A_ub=matrix,
                b_ub=[H-free(tp),free(tm)],bounds=[(-1,1)]*(len(knots)-1)+[(None,None)])
    assert res.success and res.fun>1e-7
    print("narrow slab x,v,a; two touch times; LP margin:", x,v,a,tm,tp,-res.fun)
    b,L=.2,1.5
    print("blind interval: displacement-only and necessary width:",
          b*L*L,b*L*L+2*peak(b*L,0,b)[0])


if __name__ == "__main__":
    verify()
