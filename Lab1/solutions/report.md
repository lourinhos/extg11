
## Exercise 2: Deleterious Mutations

### Why do well-known, deadly, congenital diseases persist in humans?

Our simulations show that harmful alleles never fully disappear as long as new mutations keep appearing. With a mutation rate of $u = 0.0001$ and selection against $A_2$, the allele settles at a low but non-zero frequency — a balance between mutation adding new copies and selection removing them (**mutation–selection balance**). For the semi-dominant case ($w_2 = 0.95$, $w_3 = 0.9$), selection can act on carriers too, so the equilibrium frequency is very low ($\hat{q} = 0.0019$). When the allele is fully recessive ($w_2 = 1.0$, $w_3 = 0.9$), only homozygotes are affected, so it settles at a much higher frequency ($\hat{q} = 0.031$). Even a recessive lethal ($w_3 = 0$) persists at $\hat{q} = 0.0099$. This is why diseases like cystic fibrosis persist — mutation keeps feeding in new copies, and selection can't "see" recessive alleles hiding in healthy carriers.

### Where are rare (and often deleterious) alleles often found?

Mostly in **heterozygous carriers** — people who carry one copy but are healthy. When an allele is rare, almost all copies are in heterozygotes rather than homozygotes. For example, at $q = 0.01$, only about 1 in 100 copies of $A_2$ ends up in a homozygote. So the majority of harmful alleles are carried around by people who show no symptoms.

### Why might many deleterious alleles often be strongly recessive?

Our simulations show that if a harmful allele has even a small dominant effect ($h > 0$), selection removes it quickly — the semi-dominant case settled at just $q = 0.0019$, whereas the fully recessive case reached $q = 0.031$, over 16 times higher. Recessive alleles stick around much longer because carriers are healthy and invisible to selection. Over time, the harmful alleles we still see in populations tend to be the recessive ones — the dominant ones got weeded out long ago. It's basically a filtering effect: what survives is what's hardest to remove.

### What are some implications for human health?

- **Carrier screening** matters because most disease alleles are in healthy carriers. If two carriers have a child, there's a 1/4 chance the child is affected.
- **Inbreeding** increases the chance of being homozygous, which exposes these hidden recessive alleles and raises disease risk.
- **We can't get rid of** recessive diseases through selection alone — mutation will always bring them back.
- Some harmful alleles stick around at higher-than-expected frequencies because carriers actually have a fitness benefit (e.g., sickle-cell carriers resist malaria).

---

## Exercise 3: Genetic Drift

### What are the possible outcomes of genetic drift in a finite population?

With no selection or mutation, there are only two possible endpoints: the allele either reaches a frequency of 1 (**fixation**) or drops to 0 (**loss**). Both alleles can't coexist forever — random fluctuations will eventually push one out. This is different from the infinite-population case where alleles can stay at stable intermediate frequencies.

### What is the probability that a neutral allele will go to fixation?

For a neutral allele, the fixation probability equals its starting frequency: $P(\text{fix}) = q_0$. Our simulations confirmed this — at $q_0 = 0.5$ we got $P(\text{fix}) = 0.481$, at $q_0 = 0.4$ we got $0.401$, at $q_0 = 0.1$ we got $0.105$, all close to $q_0$ itself.

### How long does it take for a neutral allele to go to fixation on average?

Starting at $q_0 = 0.5$ with $N = 100$, the average fixation time was about 274 generations. As $q_0$ decreased, fixation took longer — at $q_0 = 0.1$ it was about 395 generations. Theory predicts roughly $4N = 400$ generations for a new mutation, which lines up with what we observed.

### How long does it take for a neutral allele to go extinct on average?

Extinction is much faster than fixation. At $q_0 = 0.5$, average extinction time was 268 generations, but at $q_0 = 0.1$ it dropped to 102, and at $q_0 = 0.002$ it was just 5 generations. For $q_0 = 1/(2N)$, our simulations gave an average of about 10 generations for $N = 100$ and 27 for $N = 500$ and $N = 1000$. Most new mutations just disappear quickly by bad luck before they get a chance to spread.

### Is there a relation between these average times and N or q₀? Why?

Yes. Bigger populations mean slower drift, so everything takes longer. We saw extinction times go from ~10 generations ($N = 100$) to ~27 generations ($N = 500, 1000$) for $q_0 = 1/(2N)$. Higher $q_0$ makes fixation more likely and extinction slower. This makes sense because in a bigger population, each generation's random sampling has a smaller effect on allele frequency, so it takes more generations for frequencies to drift to the boundaries.

### How well did your simulations match analytic predictions for heterozygosity decay?

The theory predicts heterozygosity decays as $H_t = H_0(1 - \frac{1}{2N})^t$. We compared the mean $H$ at generation 500 across replicates to this prediction. With $N = 100$, the analytic value is $H_{500} = 0.041$. With 20 replicates we got 0.071 (off by 0.030), with 100 replicates we got 0.029 (off by 0.012), and with 200 replicates we got 0.044 (off by just 0.003). With $N = 500$, the analytic value is 0.303, and 100 replicates nailed it exactly (0.303). So the simulations match the theory well once you have enough replicates.

### How many replicate simulations did you need before convergence?

Around 100–200 replicates were needed. With only 20 replicates the error was noticeable (e.g. off by 0.03 for $N = 100$, off by 0.06 for $N = 500$). At 200 replicates the difference shrank to less than 0.01 in both cases.

### Did you try different population sizes? How did N affect the results?

Yes. With $N = 500$, heterozygosity at generation 500 was still around 0.30, while with $N = 100$ it had dropped to about 0.04. Larger $N$ means slower decay, exactly as the $(1 - 1/(2N))^t$ formula predicts. With larger $N$ we also needed more replicates because individual runs vary more when drift is weaker.

---

## Exercise 4: Selection in Finite Populations

### What happened as you increased population size, selection strength, or both?

The beneficial allele fixed more often. With weak selection and small $N$ ($s = 0.001$, $N = 100$), fixation probability was just 0.018 — barely above the neutral expectation of $q_0 = 0.02$. With strong selection and large $N$ ($s = 0.1$, $N = 10000$), fixation probability hit 1.0 — every single run ended in fixation.

### At what point did dynamics transition from drift-dominated to selection-dominated? What is your "rule of thumb"?

The key number is $Ns$ — population size times selection coefficient. Our rule of thumb:

- $Ns < 1$: drift wins — the allele basically behaves as if it were neutral
- $Ns > 1$: selection wins — you start to see a clear advantage for the beneficial allele

Our table shows this clearly. For $s = 0.001$, $N = 100$ ($Ns = 0.1$), $P(\text{fix}) = 0.018 \approx q_0$. For $s = 0.01$, $N = 1000$ ($Ns = 10$), $P(\text{fix}) = 0.329$, well above neutral. For $s = 0.1$, $N = 10000$ ($Ns = 1000$), $P(\text{fix}) = 1.0$.

### Do beneficial alleles always go to fixation in finite populations?

No. Even with a clear selective advantage, alleles are often lost. With $w_2 = 1.01$, $w_3 = 1.02$ and $N = 100$, only 7% of runs ended in fixation. Even with $s = 0.01$ and $N = 100$, fixation probability was just 3.5%. Most beneficial mutations get unlucky early on and disappear before selection can help them.

### What does this tell you about natural selection as a "hill-climbing" algorithm?

Evolution is not a guaranteed march toward the best outcome. Random drift can override selection — good alleles get lost and bad ones sometimes spread. Our results show that even with a 10% fitness advantage ($s = 0.1$), fixation only reaches ~32% in a population of 100. The smaller the population or the weaker the selection, the more randomness matters. Populations don't always end up at the fitness peak. Evolution is more like a wobbly random walk with a slight uphill bias, not a precise optimizer.
