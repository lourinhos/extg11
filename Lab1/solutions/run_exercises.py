import sys, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import random
from scipy import *

exec(open('./Allele-A1-functions.py').read())

print("=== EXERCISE 1 ===")

# 1a
q = determTraj(q0=0.01, w_1=0.95, w_2=1, w_3=0.95, nGen=300)
print(f"  Overdominance (w1=0.95, w2=1, w3=0.95): final q={q[-1]:.4f}")

# 1b
print("  Varying dominance with s=0.04:")
for w2 in [1.0, 1.01, 1.02, 1.03, 1.04]:
    q = determTraj(q0=0.01, w_1=1.0, w_2=w2, w_3=1.04, nGen=300)
    h = (w2 - 1) / 0.04
    print(f"    h={h:.2f} (w2={w2}): final q={q[-1]:.4f}")

# 1c
q = determTraj(q0=0.01, w_1=1.0, w_2=1.05, w_3=1.0, nGen=300)
print(f"  Het. advantage (w1=1, w2=1.05, w3=1): final q={q[-1]:.4f}")

print("\n=== EXERCISE 2 ===")
q1 = determTraj(q0=0.001, w_1=1, w_2=0.95, w_3=0.9, nGen=1000, u=0.0001)
q2 = determTraj(q0=0.001, w_1=1, w_2=1.0,  w_3=0.9, nGen=1000, u=0.0001)
q3 = determTraj(q0=0.001, w_1=1, w_2=1.0,  w_3=0.0, nGen=1000, u=0.0001)
print(f"  w2=0.95, w3=0.9 (semi-dominant): final q = {q1[-1]:.6f}")
print(f"  w2=1.0,  w3=0.9 (recessive):     final q = {q2[-1]:.6f}")
print(f"  w2=1.0,  w3=0.0 (recessive lethal): final q = {q3[-1]:.6f}")

print("\n=== EXERCISE 3 ===")
np.random.seed(42)

# 3a
print("  P(fix) for neutral allele, N=100, 1000 sims:")
for q0 in [0.5, 0.4, 0.3, 0.2, 0.1, 1/500]:
    fixations = 0
    ext_times = []
    fix_times = []
    for _ in range(1000):
        result = WFSimToEnd(q0=q0, w_1=1, w_2=1, w_3=1, N=100)
        if result[0] == 1:
            fixations += 1
            fix_times.append(result[1])
        else:
            ext_times.append(result[1])
    avg_ext = np.mean(ext_times) if ext_times else float('nan')
    avg_fix = np.mean(fix_times) if fix_times else float('nan')
    print(f"    q0={q0:.4f}: P(fix)={fixations/1000:.3f}, avg_t_ext={avg_ext:.0f}, avg_t_fix={avg_fix:.0f}")


#3b
print("\n  Avg time to extinction for q0=1/(2N):")
for N in [100, 500, 1000]:
    q0 = 1/(2*N)
    times = []
    for _ in range(1000):
        result = WFSimToEnd(q0=q0, w_1=1, w_2=1, w_3=1, N=N)
        times.append(result[1])
    print(f"    N={N}, q0={q0:.5f}: avg_time={np.mean(times):.1f}")

# =============================================================================
# EXERCISE 4: Selection in finite populations
# =============================================================================
print("\n=== EXERCISE 4 ===")

# 4a. Initial test: q0=0.02, w2=1.01, w3=1.02, 1000 sims
fixations = 0
for _ in range(1000):
    result = WFSimToEnd(q0=0.02, w_1=1.0, w_2=1.01, w_3=1.02, N=100)
    fixations += result[0]
print(f"  q0=0.02, w2=1.01, w3=1.02, N=100, 1000 sims: P(fix)={fixations/1000:.3f}")

# 4b. Full table: q0=0.02, 1000 sims per cell (10000 is very slow for N=10000)
nSims = 1000
print(f"\n  Fixation table (q0=0.02, {nSims} sims per cell):")
print(f"  {'Fitnesses':<30} {'N=100':>8} {'N=1000':>8} {'N=10000':>8}")
fitness_sets = [
    (1.0005, 1.001),
    (1.005,  1.01),
    (1.05,   1.1),
]
for w2, w3 in fitness_sets:
    row = f"  w2={w2}, w3={w3}"
    row = f"{row:<30}"
    for N in [100, 1000, 10000]:
        fixations = 0
        for _ in range(nSims):
            result = WFSimToEnd(q0=0.02, w_1=1.0, w_2=w2, w_3=w3, N=N)
            fixations += result[0]
        row += f" {fixations/nSims:.4f} "
    print(row)

print("\nDone.")
