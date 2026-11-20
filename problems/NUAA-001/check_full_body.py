"""Numerical checks in unscaled 3-D quaternion, velocity and four-motor coordinates."""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.spatial.transform import Rotation


def hat(x):
    return np.cross(x[:, None], np.eye(3), axisa=0, axisb=0, axisc=0)


mixer = np.array([[1, 1, 1, 1], [1, -1, -1, 1],
                  [-1, -1, 1, 1], [0.7, -0.7, 0.7, -0.7]])
inertia = np.array([0.8, 1.0, 1.2])
decay = np.array([1.1, 1.25, 1.4, 1.5])
hover = np.full(4, 0.25)
initial = np.r_[np.zeros(6), [0, 0, 0, 1], np.zeros(3), hover]
e3 = np.array([0., 0., 1.])


def dynamics(state, commands):
    quat, omega, motors = state[6:10], state[10:13], state[13:]
    rotor = Rotation.from_quat(quat)
    thrust, *moments = mixer @ motors
    qdot = np.r_[quat[3]*omega + np.cross(quat[:3], omega), -quat[:3]@omega]/2
    return np.r_[state[3:6], thrust*rotor.apply(e3)-e3, qdot,
                 (moments-np.cross(omega, inertia*omega))/inertia,
                 decay*(commands-motors)]


def endpoint(h, command, start=initial):
    sol = solve_ivp(lambda t, x: dynamics(x, command(t/h)), (0, h), start,
                    method="DOP853", rtol=1e-12, atol=1e-15)
    assert sol.success
    return sol.y[:, -1]


def scaled(state, h):
    return np.r_[state[:2]/h**5, state[2]/h**3,
                 state[3:5]/h**4, state[5]/h**2,
                 Rotation.from_quat(state[6:10]).as_rotvec()/h**3,
                 state[10:13]/h**2, (state[13:]-hover)/h]


def short_window():
    control = np.array([.08, -.11, .12, -.05])
    torque_rate = (mixer[1:]@(decay*control))/inertia
    prediction = np.cross(torque_rate, e3)[:2]/120
    errors = []
    for h in (.2, .1, .05, .025):
        end = endpoint(h, lambda s: hover+control)
        errors.append(np.linalg.norm(end[:2]/h**5-prediction))
    assert all(a*0.57 > b for a, b in zip(errors, errors[1:])), errors
    # Independent endpoint Jacobian by raw nonlinear integration. First two
    # torque channels need 5 moments, the vertical and yaw channels need 3.
    inv = np.linalg.inv(mixer@np.diag(decay))
    results = []
    for h in (.3, .1):
        cols = []
        for channel, dimension in enumerate((3, 5, 5, 3)):
            for degree in range(dimension):
                perturb = inv[:, channel]*.03
                plus = endpoint(h, lambda s: hover+perturb*(1-s)**degree)
                minus = endpoint(h, lambda s: hover-perturb*(1-s)**degree)
                cols.append((scaled(plus, h)-scaled(minus, h))/.06)
        jac = np.array(cols).T
        singular_values = np.linalg.svd(jac, compute_uv=False)
        assert len(singular_values) == 16 and singular_values[-1] > 1e-8
        results.append(singular_values[-1])
    print("3D h^5 support errors:", errors, "scaled endpoint sigma_min:", results)


def hover_linearization(thrust):
    A = np.zeros((16, 16))
    A[:3, 3:6] = A[6:9, 9:12] = np.eye(3)
    A[3:6, 6:9] = -thrust*hat(e3)
    A[3:6, 12:] = np.outer(e3, mixer[0])
    A[9:12, 12:] = mixer[1:]/inertia[:, None]
    A[12:, 12:] = -np.diag(decay)
    B = np.zeros((16, 4))
    B[12:] = np.diag(decay)
    return A, B


def recovery():
    A, B = hover_linearization(1)
    # Only middle half of the window is allowed to correct. Held motor
    # commands are integrated exactly in the model used to design the gain.
    dt, slots = 0.125, 24
    H = dt*slots
    augmented = np.zeros((20, 20))
    augmented[:16, :16], augmented[:16, 16:] = A, B
    transition = expm(augmented*dt)
    F, U = transition[:16, :16], transition[:16, 16:]
    active = list(range(6, 18))
    effect = np.hstack([np.linalg.matrix_power(F, slots-1-j)@U for j in active])
    assert np.linalg.matrix_rank(effect) == 16
    zero, B0 = hover_linearization(0)
    reachable = np.hstack([np.linalg.matrix_power(zero, j)@B0 for j in range(16)])
    assert np.linalg.matrix_rank(reachable) == 12
    initial_error = np.array([1, -1, .25]+[0]*13, dtype=float)
    errors = []
    for size in (1e-4, 5e-5, 2.5e-5):
        rhs = -expm(A*H)@(size*initial_error)
        corrections = np.linalg.lstsq(effect, rhs, rcond=None)[0].reshape(-1, 4)
        assert np.max(abs(corrections)) < .2
        state = initial.copy()
        state[:3] += size*initial_error[:3]
        for j in range(slots):
            desired = hover + (corrections[j-active[0]] if j in active else 0)
            state = endpoint(dt, lambda s: desired, state)
        errors.append(np.linalg.norm(scaled(state, 1)))
    assert all(b < .28*a for a, b in zip(errors, errors[1:])), errors
    print("held-command correction, raw nonlinear endpoint errors:", errors)


if __name__ == "__main__":
    short_window()
    recovery()
