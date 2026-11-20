"""Compare matched rotor commands in a full rigid-body simulation."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy.spatial.transform import Rotation
from check_full_body import decay, dynamics, hover, initial


def switches(z):
    # Equivalent basis 1,t,(exp(zt)-1-zt)/z², evaluated by quadrature.
    abscissa, weights = np.polynomial.legendre.leggauss(10)

    def residual(cuts):
        result = np.zeros(3)
        edges = np.r_[0, cuts, 1]
        for j, (a, b) in enumerate(zip(edges, edges[1:])):
            t = (a+b)/2 + (b-a)*abscissa/2
            basis = np.array([0*t+1, t, (np.expm1(z*t)-z*t)/z**2])
            result += (-1)**j*(b-a)/2*(basis @ weights)
        return result

    sol = root(residual, [.15, .5, .85], tol=1e-11)
    assert np.max(abs(residual(sol.x))) < 1e-12
    assert np.all(np.diff(np.r_[0, sol.x, 1]) > 0)
    return sol.x


def compare(n):
    local = np.array([switches(rate/n) for rate in decay])
    y = initial.copy()
    y[6:10] = Rotation.from_rotvec([.25, -.1, .2]).as_quat()
    y[10:13] = [.2, .13, -.06]
    pair = np.r_[y, y]
    errors = np.zeros(3)
    for k in range(n):
        times = np.unique(np.r_[0, local.flatten(), 1])
        for a, b in zip(times, times[1:]):
            parity = np.sum(local < (a+b)/2, axis=1)
            alternate = hover + .25*(-1.)**(parity + np.array([0, 1, 1, 0]))

            def rhs(t, states):
                vals = []
                for state, command in zip(states.reshape(2, -1), (hover, alternate)):
                    derivative = dynamics(state, command)
                    derivative[3] += .1*(-1)**(k//(n//4))
                    derivative[11] += .05*(-1)**(k//(n//2))
                    vals.append(derivative)
                return np.concatenate(vals)

            start, stop = (k+a)/n, (k+b)/n
            sol = solve_ivp(rhs, (start, stop), pair,
                            t_eval=np.linspace(start, stop, 5), method="DOP853",
                            rtol=1e-11, atol=1e-14)
            assert sol.success
            pair = sol.y[:, -1]
            ref, changed = sol.y.reshape(2, 17, -1)
            position = np.linalg.norm(ref[:3]-changed[:3], axis=0).max()
            left, right = (Rotation.from_quat(block[6:10].T) for block in (ref, changed))
            angle = (left.inv()*right).magnitude().max()
            errors[:2] = np.maximum(errors[:2], [position, angle])
        errors[2] = max(errors[2], np.max(abs(pair[13:17]-pair[30:34])))
    return errors


if __name__ == "__main__":
    values = np.array([compare(n) for n in (4, 8, 16)])
    print("columns: position, rotation angle, motor grid error")
    print(values)
    assert np.max(values[:, 2]) < 1e-11
    assert np.all(values[1:, :2] < .17*values[:-1, :2])
