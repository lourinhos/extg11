import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define the interspecific competition model
def icomp(t, n):
    global r, k, alfa, beta
    dn1_dt = r[0] * n[0] * ((k[0] - n[0] - alfa * n[1]) / k[0])
    dn2_dt = r[1] * n[1] * ((k[1] - n[1] - beta * n[0]) / k[1])
    return [dn1_dt, dn2_dt]

def phase_diagram(dat):
    global r, k, alfa, beta
    plt.figure(2)
    plt.plot(dat[0, :], dat[1, :], label='Population sizes')
    plt.xlabel('Species 1')
    plt.ylabel('Species 2')
    plt.plot(dat[0, 0], dat[1, 0], '*', label='Starting point', clip_on=False, zorder=3)  # Mark the starting point

    # Isoclines
    y1 = np.linspace(0, k[0] / alfa, 10)
    x1 = k[0] - alfa * y1
    x2 = np.linspace(0, k[1] / beta, 10)
    y2 = k[1] - beta * x2
    plt.plot(x1, y1, 'r', label='dn1/dt=0')
    plt.plot(x2, y2, 'g', label='dn2/dt=0')
    plt.legend()
    plt.axis([0, max(max(x1), max(x2)), 0, max(max(y1), max(y2))])
    plt.title('Phase Diagram')


def run_case(case_name, params, n0, t0, tend, file_prefix):
    global r, k, alfa, beta
    r = params['r']
    k = params['k']
    alfa = params['alfa']
    beta = params['beta']

    t_eval = np.linspace(t0, tend, 1000)
    sol = solve_ivp(icomp, [t0, tend], n0, t_eval=t_eval)

    t = sol.t
    dat = sol.y

    plt.figure(figsize=(7, 4))
    plt.plot(t, dat[0, :], label='Species 1')
    plt.plot(t, dat[1, :], label='Species 2')
    plt.plot(t, np.zeros_like(t), 'k--', label='Zero Line')
    plt.xlabel('Time')
    plt.ylabel('Population Density')
    plt.legend()
    plt.title(case_name)
    plt.tight_layout()
    plt.savefig(f'{file_prefix}_time.png', dpi=200)
    plt.close()

    plt.figure(figsize=(6, 5))
    plt.plot(dat[0, :], dat[1, :], label='Population sizes')
    plt.xlabel('Species 1')
    plt.ylabel('Species 2')
    plt.plot(dat[0, 0], dat[1, 0], '*', label='Starting point', clip_on=False, zorder=3)

    y1 = np.linspace(0, k[0] / alfa, 10)
    x1 = k[0] - alfa * y1
    x2 = np.linspace(0, k[1] / beta, 10)
    y2 = k[1] - beta * x2
    plt.plot(x1, y1, 'r', label='dn1/dt=0')
    plt.plot(x2, y2, 'g', label='dn2/dt=0')
    plt.legend()
    plt.axis([0, max(max(x1), max(x2)), 0, max(max(y1), max(y2))])
    plt.title(f'{case_name} Phase Diagram')
    plt.tight_layout()
    plt.savefig(f'{file_prefix}_phase.png', dpi=200)
    plt.close()


if __name__ == "__main__":
    t0, tend = 0, 300

    cases = [
        (
            '28 Coexistence',
            {
                'r': np.array([1.0, 1.0]),
                'k': np.array([1000.0, 1000.0]),
                'alfa': 0.5,
                'beta': 0.5,
            },
            [1, 1],
            'case28',
        ),
        (
            '29 Species 2 Wins',
            {
                'r': np.array([1.0, 1.0]),
                'k': np.array([1000.0, 1000.0]),
                'alfa': 1.5,
                'beta': 0.5,
            },
            [1, 1],
            'case29',
        ),
        (
            '30 Species 1 Wins',
            {
                'r': np.array([1.0, 1.0]),
                'k': np.array([1000.0, 1000.0]),
                'alfa': 0.5,
                'beta': 1.5,
            },
            [1, 1],
            'case30',
        ),
        (
            '31 Outcome A',
            {
                'r': np.array([1.0, 1.0]),
                'k': np.array([1000.0, 1000.0]),
                'alfa': 1.5,
                'beta': 1.5,
            },
            [800, 100],
            'case31_a',
        ),
        (
            '31 Outcome B',
            {
                'r': np.array([1.0, 1.0]),
                'k': np.array([1000.0, 1000.0]),
                'alfa': 1.5,
                'beta': 1.5,
            },
            [100, 800],
            'case31_b',
        ),
    ]

    for case_name, params, n0, file_prefix in cases:
        run_case(case_name, params, n0, t0, tend, file_prefix)

    plt.show()

