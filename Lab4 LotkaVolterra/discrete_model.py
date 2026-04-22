from matplotlib import pyplot as plt
import numpy as np
from population_growth import intra_comp

def discrete_intra_comp(n0, r, k, timesteps):
    n = np.zeros(timesteps + 1)
    n[0] = n0

    for t in range(timesteps):
        n[t + 1] = n[t] + (r * n[t] * (1 - n[t] / k))  # Equation 16

    t_values = np.arange(0, timesteps + 1)

    nt_cont = intra_comp(n0, r, t_values, k, plot=False)  # Continuous model for comparison

    plt.figure()
    plt.plot(t_values, n, label='Discrete Model')
    plt.plot(t_values, nt_cont, label='Continuous Model')
    plt.xlabel('Time Steps')
    plt.ylabel('Population Size')
    plt.title('Discrete Time Model')
    plt.legend()
    plt.savefig('discrete_intra_comp.png')  # Save the figure as a PNG file
    plt.show()

if __name__ == "__main__":
    n0 = 200  # Initial population size
    r = 1.2 # Growth rate
    k = 1000  # Carrying capacity

    timesteps = 100  # Number of time steps

    discrete_intra_comp(n0, r, k, timesteps)
    