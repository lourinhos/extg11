from matplotlib import pyplot as plt
import numpy as np

def growth_model(n0, r, t):
    nt = n0 * np.exp(r * t)  # Growth equation

    plt.figure()
    plt.plot(t, nt, color='blue', label=f'r = {r}')
    plt.legend()
    plt.xlabel('Time (years)')
    plt.ylabel('Population Size')
    plt.title('Continuous Population Growth')
    plt.show()


def intra_comp(n0, r, t, k, plot=True):
    nt = k / (1 + (k / n0 - 1) * np.exp(-r * t))  

    if plot:
        plt.figure()
        plt.plot(t, nt)
        plt.xlabel('Time (years)')
        plt.ylabel('Population Size')
        plt.title('Intraspecific Competition')
        plt.savefig('intra_comp0.png')  # Save the figure as a PNG file
        plt.show()

    return nt

def intra_comp_growth(n, r, k):
    dndt = r * n * (1 - n / k)  

    plt.figure()
    plt.plot(n, dndt)
    plt.xlabel('Population Size')
    plt.ylabel('dN/dt')
    plt.title('Intraspecific Competition Growth')
    plt.savefig('intra_comp_growth.png')  # Save the figure as a PNG file
    plt.show()



if __name__ == "__main__":
    n = np.arange(0, 1010, 10) # Population size range
    r = 1 # Growth rate
    k = 1000 # Carrying capacity
    intra_comp_growth(n, r, k)

