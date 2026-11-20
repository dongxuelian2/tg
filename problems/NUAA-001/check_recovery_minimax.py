"""Paired nonlinear replays at the same torque, with two distinct force histories."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation
from check_full_body import dynamics, initial, hover


def replay(windy):
    left = initial.copy()
    right = initial.copy()
    left[6:10] = Rotation.from_rotvec([.12, -.18, .03]).as_quat()
    right[6:10] = Rotation.from_rotvec([.13, -.20, .05]).as_quat()
    right[:3] = [.03, -.02, .01]
    state = np.r_[left, right]
    differences = []
    for k in range(10):
        gust = np.array([3.0*(-1)**k, -2.1*(k%3-1), 1.5*np.cos(k)]) if windy else np.zeros(3)
        torque_acc = np.array([.02*np.cos(k), -.04, .03*(-1)**k])
        def rhs(t, y):
            out = []
            for j, body in enumerate(y.reshape(2, -1)):
                delta = (.012*np.array([1., -1., 1., -1.]) if j else np.zeros(4))
                cmd = hover + .04*np.sin(t+np.arange(4)) + delta
                val = dynamics(body, cmd)
                val[3:6] += gust
                val[10:13] += torque_acc
                out.append(val)
            return np.concatenate(out)
        times = np.linspace(k/10, (k+1)/10, 7)
        sol = solve_ivp(rhs, (times[0], times[-1]), state, t_eval=times,
                        rtol=1e-12, atol=1e-14, method="DOP853")
        assert sol.success
        left, right = sol.y.reshape(2, 17, -1)
        differences.extend((right-left).T)
        state = sol.y[:, -1]
    return np.array(differences), state


if __name__ == "__main__":
    calm, yc = replay(False)
    windy, yw = replay(True)
    error = np.max(abs(calm-windy))
    actual = np.linalg.norm(yc[:3]-yw[:3])
    assert error < 1e-10 and actual > .1
    print("absolute position change / change in pairwise trajectory difference:", actual, error)
