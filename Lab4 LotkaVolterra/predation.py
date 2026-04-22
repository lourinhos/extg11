import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def pred(t, n):
    global r, a, b, k, d
    n1, n2 = max(n[0], 0), max(n[1], 0)  # Ensure non-negative populations
    dn1_dt = r * n1 * (1 - n1 / k) - a * n1 * n2
    dn2_dt = b * a * n1 * n2 - d * n2
    return [dn1_dt, dn2_dt]

if __name__ == "__main__":
    # Initial conditions
    n0 = [10, 10]  # Density of prey and predator at t=0
    t0, tend = 0, 400  # Simulation start and stop
    t_eval = np.linspace(t0, tend, 1000)

    cases = [
        (
            'a) Prey goes extinct',
            {'r': 0.2, 'a': 0.15, 'b': 0.08, 'd': 0.02, 'k': 2000},
        ),
        (
            'b) Predator goes extinct',
            {'r': 0.5, 'a': 0.005, 'b': 0.01, 'd': 0.5, 'k': 2000},
        ),
        (
            'c) Damped oscillations to coexistence',
            {'r': 0.6, 'a': 0.01, 'b': 0.02, 'd': 0.1, 'k': 2000},
        ),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(18, 9))

    for col, (case_title, params) in enumerate(cases):
        global r, a, b, k, d
        r = params['r']
        a = params['a']
        b = params['b']
        d = params['d']
        k = params['k']

        sol = solve_ivp(pred, [t0, tend], n0, t_eval=t_eval)
        t = sol.t
        dat = sol.y

        density_ax = axes[0, col]
        density_ax.plot(t, dat[0, :], label='Prey')
        density_ax.plot(t, dat[1, :], label='Predator')
        density_ax.set_xlabel('Time')
        density_ax.set_ylabel('Population Density')
        density_ax.set_title(
            f"{case_title}\nr={params['r']}, a={params['a']}, b={params['b']}, d={params['d']}, k={params['k']}"
        )
        density_ax.legend()

        phase_ax = axes[1, col]
        phase_ax.plot(dat[0, :], dat[1, :])
        phase_ax.set_xlabel('Prey Density')
        phase_ax.set_ylabel('Predator Density')
        phase_ax.set_title('Phase Diagram')

    fig.suptitle('Predator-Prey Dynamics', fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    plt.show()


