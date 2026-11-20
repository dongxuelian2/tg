"""Integrate both drag models and recover endpoint data only from relative telemetry."""
import numpy as np
from scipy.integrate import solve_ivp


lags = np.array([.7, 1.3, 1.8])
initial = np.arange(18.).reshape(3, 6)/30


def check(air, freq):
    def rhs(t, flat):
        y = flat.reshape(3, 9)
        accel = np.cos(freq*t + np.arange(9).reshape(3, 3))
        wind = np.sin(t*np.array([1., 2., 3.]))
        ground = accel - lags[:, None]*y[:, 3:6]
        ground += lags[:, None]*wind if air else wind
        return np.hstack([y[:, 3:6], ground, accel]).ravel()

    start = np.c_[initial, np.zeros((3, 3))]
    sol = solve_ivp(rhs, (0, 1.7), start.ravel(), rtol=1e-11, atol=1e-13,
                    method="DOP853")
    assert sol.success
    delta = sol.y[:, -1].reshape(3, 9)-start
    rel = delta[1:]-delta[0]
    blocks = rel[:, 6:]-rel[:, 3:6]-lags[1:, None]*rel[:, :3]
    recovered = np.linalg.lstsq((lags[1:]-lags[0])[:, None], blocks, rcond=None)[0][0]
    if air:
        recovered = delta[0, 6:]-lags[0]*recovered
    expected = delta[0, 3:6] if air else delta[0, :3]
    assert np.max(abs(recovered-expected)) < 1e-9


if __name__ == "__main__":
    for air in [False, True]:
        for freq in [.3, 1.7, 5.]:
            check(air, freq)
    mat = np.array([[np.log(2)]*3, [.5,.25,.125], [.375,.09375,.0234375]])
    coeff = np.linalg.solve(mat, np.ones(3))
    coeff /= np.max(abs(coeff))
    impulse, v0, v1 = mat @ coeff
    assert abs(v0-v1) < 1e-12 and abs(impulse-v0) < 1e-12 and abs(v0) > .01
    print("Six integral reconstructions passed; indistinguishable terminal speed:", v0)
