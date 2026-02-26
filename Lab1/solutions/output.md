### EXERCISE 1 
  Overdominance (w1=0.95, w2=1, w3=0.95): final q=0.4989
  Varying dominance with s=0.04:
    h=0.00 (w2=1.0): final q=0.0113
    h=0.25 (w2=1.01): final q=0.2272
    h=0.50 (w2=1.02): final q=0.7851
    h=0.75 (w2=1.03): final q=0.8675
    h=1.00 (w2=1.04): final q=0.8435
  Het. advantage (w1=1, w2=1.05, w3=1): final q=0.4984

### EXERCISE 2 
  w2=0.95, w3=0.9 (semi-dominant): final q = 0.001900
  w2=1.0,  w3=0.9 (recessive):     final q = 0.031404
  w2=1.0,  w3=0.0 (recessive lethal): final q = 0.009901

### EXERCISE 3
  P(fix) for neutral allele, N=100, 1000 sims:
    q0=0.5000: P(fix)=0.481, avg_t_ext=268, avg_t_fix=274
    q0=0.4000: P(fix)=0.401, avg_t_ext=231, avg_t_fix=296
    q0=0.3000: P(fix)=0.314, avg_t_ext=219, avg_t_fix=335
    q0=0.2000: P(fix)=0.198, avg_t_ext=161, avg_t_fix=342
    q0=0.1000: P(fix)=0.105, avg_t_ext=102, avg_t_fix=395
    q0=0.0020: P(fix)=0.001, avg_t_ext=5, avg_t_fix=231

  Avg time to extinction for q0=1/(2N):
    N=100, q0=0.00500: avg_time=10.2
    N=500, q0=0.00100: avg_time=27.0
    N=1000, q0=0.00050: avg_time=27.4

### EXERCISE 4
  q0=0.02, w2=1.01, w3=1.02, N=100, 1000 sims: P(fix)=0.070

  Fixation table (q0=0.02, 1000 sims per cell):
  Fitnesses                         N=100   N=1000  N=10000
  w2=1.0005, w3=1.001          0.0180  0.0530  0.3220 
  w2=1.005, w3=1.01            0.0350  0.3290  0.9830 
  w2=1.05, w3=1.1              0.3210  0.9790  1.0000 

