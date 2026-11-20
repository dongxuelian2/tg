"""Compare two-segment predictions against the unreset rigid-body trajectory."""
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp
from check_full_body import dynamics, hover_linearization, hover, initial, scaled


def run():
    dt, slots = 0.125, 24
    A, B = hover_linearization(1)
    # One shared scalar lateral wind. No disturbance torque is available.
    augmented = np.zeros((21, 21))
    augmented[:16, :16] = A
    augmented[:16, 16:20] = B
    augmented[3, 20] = 1
    step = expm(dt * augmented)
    F, U, wind = step[:16, :16], step[:16, 16:20], step[:16, 20]
    active = range(6, 18)
    effect = np.hstack([np.linalg.matrix_power(F, slots-1-j) @ U for j in active])
    gain = np.linalg.lstsq(effect, -np.linalg.matrix_power(F, slots), rcond=None)[0]
    assert np.linalg.norm(effect @ gain + np.linalg.matrix_power(F, slots)) < 1e-9
    errors = []
    for amplitude in (2e-5, 1e-5):
        state = initial.copy()
        previous_wind_endpoint = np.zeros(16)
        worst = 0
        for index, sign in enumerate((1, -1, -1, 1, -1, 1)):
            actual_inputs = (gain @ scaled(state, 1)).reshape(-1, 4)
            predicted_inputs = (gain @ previous_wind_endpoint).reshape(-1, 4)
            predicted = previous_wind_endpoint.copy()
            current_wind = np.zeros(16)
            for j in range(slots):
                gust = amplitude * sign * (1 if j < 9 + index else -0.4)
                command = hover + (actual_inputs[j-active.start] if j in active else 0)
                assert np.all((command > 0) & (command < 1))

                def rhs(t, value):
                    derivative = dynamics(value, command)
                    derivative[3] += gust
                    return derivative

                sol = solve_ivp(rhs, (0, dt), state, method="DOP853",
                                rtol=1e-12, atol=1e-14)
                assert sol.success
                state = sol.y[:, -1]
                predicted = F @ predicted + wind * gust
                if j in active:
                    predicted += U @ predicted_inputs[j-active.start]
                current_wind = F @ current_wind + wind * gust
                worst = max(worst, np.linalg.norm(scaled(state, 1) - predicted))
            previous_wind_endpoint = current_wind
        errors.append(worst)
    assert 0.23 < errors[1]/errors[0] < 0.27, errors
    print("six segments, maximum 16-state error:", errors)


if __name__ == "__main__":
    run()
