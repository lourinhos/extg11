
## Exercise 2: Deleterious Mutations

### Why do well-known, deadly, congenital diseases persist in humans?

Our simulations show that harmful alleles never fully disappear as long as new mutations keep appearing. With a mutation rate of $u = 0.0001$ and selection against $A_2$, the allele settles at a low but non-zero frequency — a balance between mutation adding new copies and selection removing them (**mutation–selection balance**). When the harmful allele is semi-dominant ($w_2 = 0.95$, $w_3 = 0.9$), selection can act on carriers too, so the equilibrium frequency is very low (~0.001). When the allele is fully recessive ($w_2 = 1.0$, $w_3 = 0.9$), only homozygotes are affected, so it settles at a higher frequency (~0.033). Even a recessive lethal ($w_3 = 0$) hangs around at $\hat{q} \approx 0.01$. This is why diseases like cystic fibrosis persist — mutation keeps feeding in new copies, and selection can't "see" recessive alleles hiding in healthy carriers.

### Where are rare (and often deleterious) alleles often found?

Mostly in **heterozygous carriers**, people who carry one copy but are healthy. When an allele is rare, almost all copies are in heterozygotes rather than homozygotes. For example, at $q = 0.01$, only about 1 in 100 copies of $A_2$ ends up in a homozygote. So the majority of harmful alleles are carried around by people who show no symptoms.

### Why might many deleterious alleles often be strongly recessive?

Our simulations show that if a harmful allele has even a small dominant effect ($h > 0$), selection removes it quickly. Recessive alleles ($h = 0$) stick around much longer because carriers are healthy and invisible to selection. Over time, the harmful alleles that we still see in populations tend to be the recessive ones — the dominant harmful ones got weeded out long ago. It's basically a filtering effect: what survives is what's hardest to remove.

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

For a neutral allele, the fixation probability simply equals its starting frequency: $P(\text{fix}) = q_0$. Our simulations matched this well — starting at $q_0 = 0.5$, about half the runs fixed; at $q_0 = 0.1$, about 10% fixed, and so on.

### How long does it take for a neutral allele to go to fixation on average?

For a brand new mutation ($q_0 = 1/(2N)$), the average time to fixation is roughly $4N$ generations. So in a population of $N = 100$, that's about 400 generations. Starting at $q_0 = 0.5$, fixation takes a few hundred generations with $N = 100$, which is what we saw in our simulations.

### How long does it take for a neutral allele to go extinct on average?

Extinction is much faster than fixation. A new mutation ($q_0 = 1/(2N)$) typically goes extinct in roughly $2\ln(2N)$ generations — that's only about 10–14 generations for $N = 100$ to $1000$. Most new mutations just disappear quickly by bad luck before they get a chance to spread.

### Is there a relation between these average times and N or q₀? Why?

Yes. Bigger populations mean slower drift, so everything takes longer — both fixation and extinction times go up with $N$. This makes sense because in a bigger population, each generation's random sampling has a smaller effect on allele frequency. Higher starting frequency $q_0$ makes fixation more likely (and faster, relatively), while low $q_0$ makes extinction quick. Our simulations showed this clearly when we compared $N = 100$ vs. $N = 1000$.

### How well did your simulations match analytic predictions for heterozygosity decay?

The theory predicts heterozygosity decays as $H_t = H_0(1 - \frac{1}{2N})^t$. When we averaged over many simulation replicates, the mean heterozygosity tracked this curve nicely. With too few replicates the average was jumpy, but it smoothed out as we added more.

### How many replicate simulations did you need before convergence?

Around 100–200 replicates gave a good match to the theoretical curve for $N = 100$. With only 10–20 replicates, the average was still quite noisy.

### Did you try different population sizes? How did N affect the results?

Yes. Larger $N$ means heterozygosity decays more slowly, which matches the formula — the $(1 - 1/(2N))^t$ term shrinks more slowly when $N$ is big. With larger $N$, we also needed more replicates to get a smooth average, since individual runs take longer to reach fixation or loss.

---

## Exercise 4: Selection in Finite Populations

### What happened as you increased population size, selection strength, or both?

As we made $N$ or $s$ bigger, the beneficial allele fixed more often. With weak selection in small populations (e.g., $s = 0.001$, $N = 100$), the fixation rate was barely above what you'd expect for a neutral allele. With strong selection in large populations (e.g., $s = 0.1$, $N = 10000$), fixation was much more common and the allele frequency curves looked smooth and predictable, like in Exercise 1.

### At what point did dynamics transition from drift-dominated to selection-dominated? What is your "rule of thumb"?

The key number is $Ns$ — population size times selection coefficient. Our rule of thumb:

- $Ns < 1$: drift wins — the allele basically behaves as if it were neutral
- $Ns > 1$: selection wins — you start to see a clear advantage for the beneficial allele

For example, $s = 0.01$ with $N = 100$ gives $Ns = 1$, which is right at the boundary. Bump $N$ up to 10000 ($Ns = 100$) and selection clearly dominates.

### Do beneficial alleles always go to fixation in finite populations?

No. Beneficial alleles are lost all the time, especially when they're rare. Even a mutation with a 5% advantage only fixes about 10% of the time when it starts as a single copy. Most beneficial mutations just get unlucky early on and disappear before selection can really help them.

### What does this tell you about natural selection as a "hill-climbing" algorithm?

It tells us that evolution is not a guaranteed march toward the best possible outcome. Random drift can override selection — good alleles get lost and bad ones sometimes spread. The smaller the population or the weaker the selection, the more randomness matters. So populations don't always end up at the fitness peak. Evolution is more like a wobbly random walk with a slight uphill bias, not a precise optimizer.
