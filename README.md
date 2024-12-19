# Fast Radio Burst Cosmology
###### Under the supervision of [Professor Victoria Kaspi](https://www.physics.mcgill.ca/~vkaspi/) at [McGill University](https://www.mcgill.ca/), and in collaboration with [Dr. Aaron Pearlman](https://sites.astro.caltech.edu/~apearlman/). In collaboration with the [CHIME experiment](https://chime-experiment.ca/en) and [Canada Compute - Cedar](https://docs.alliancecan.ca/wiki/Cedar).

![alt text](https://github.com/IsolatedSingularity/FRB-Analysis/blob/main/Plots/ACF_Sum%3B%20Fit.png)

## Objective

In quantum gravity it is possible to describe a black hole quantum tunneling into a while hole in the classical limit [[1]](https://inspirehep.net/literature/1403709). This could occur in events where a black hole collapses faster than the Hawking radiation time, and so the matter in the black hole interior must tunnel to corresponding white hole interior. These tunnelling events are thought to be connected to with periodic radio signals occuring at low time scales, known as **fast radio bursts** (FRB) [[2]](https://www.researchgate.net/publication/265730144_Fast_Radio_Bursts_and_White_Hole_Signals). These signals are usually seen to be emitted by pulsars or merger systems and as the signal travels through the interstellar medium their signals experience dispersion causing different frequencies to arrive at different times. Moreover each are subject to spacetime dependent background noise given by primordial perturbations in the very early universe. This makes it interesting to study correlations of signals at different spatial separations as it give us information on how two sources of a signal can interact with one another. This is determined by statistical osbverables such as **decorrelation and scintillation bandwidths** which respectively quantify the range of frequencies over which the phases of a radio wave emitted by a source lose correlation during propagation and the rapid fluctuations in the amplitude of the radio signal during so.

We characterize this information for fast radio burst (FRB) repeater signals detected from the direction of the M81 galaxy group. Furthermore, we numerically computed scintillation and decorrelation bandwidths of signals with the CHIME/Pulsar and CHIME/FRB systems amongst non-linear cosmological noise after performing functional time series. We also analytically looked at modelling these events with with black-white hole tunneling amplitudes for corrections to the underlying spacetime geometry. 


## Code Functionality

*readScript_.py* characterizes interaction between the fast radio bursts their their interstellar medium and their potential connection to quantum gravity phenomena. The core of the code lies in its auto-correlated functionals (ACFs) of frequencies at different points in spacetime. It offers multiple methods to compute ACFs, each tailored to different analysis requirements. These ACFs provide insights into how FRB signals evolve as they propagate through spacetime. To enhance data quality, it identifies and removes channels with amplitudes exceeding defined thresholds. The code thereafter estimates the decorrelation bandwidth, a key parameter characterizing the dynamic behavior of FRB signals in different frequency bands. It generates numerous plots to visualize the results in the form of interacting signal densities and how the scintillation bandwidth can be computed from the overlapping densities.

## Caveats

One caveat to consider is the potential behavior of the frequency functional at extreme limits. In cases where the functional approaches singularities, the code has degeneracies in the information bordering the signal distribution. Removing these degeneracies would avoid redudnant computation time on the cluster which would essential for scalability/time complexity.

## Next Steps

At this point, the problem statement has been defined and the foundational code base has been laid. Potential additions to the code would include higher dimensional black-white wormhole tunneling events with sources that are causally disconnected. It would also be informative to study this in the context of a foliated spacetime, such as in the ADM formalism.

# Fast Radio Burst Cosmology
###### Under the supervision of [Professor Victoria Kaspi](https://www.physics.mcgill.ca/~vkaspi/) at [McGill University](https://www.mcgill.ca/), and in collaboration with [Dr. Aaron Pearlman](https://sites.astro.caltech.edu/~apearlman/). In collaboration with the [CHIME experiment](https://chime-experiment.ca/en) and [Canada Compute - Cedar](https://docs.alliancecan.ca/wiki/Cedar).

![Autocorrelation Function Fit](https://github.com/IsolatedSingularity/FRB-Analysis/blob/main/Plots/ACF_Sum%3B%20Fit.png?raw=true)

## Objective

Fast Radio Bursts (FRBs) are highly energetic radio signals that last only milliseconds, with origins spanning compact objects like pulsars or binary mergers. Recently, FRBs have been theorized to connect to quantum gravity phenomena, specifically the quantum tunneling of black holes into white holes [[1]](https://inspirehep.net/literature/1403709). This tunneling could manifest as periodic bursts of radio waves due to the black hole's collapse faster than its Hawking radiation time [[2]](https://www.researchgate.net/publication/265730144_Fast_Radio_Bursts_and_White_Hole_Signals).

FRBs travel through the interstellar medium (ISM), experiencing dispersion (frequency-dependent arrival times) and scintillation (rapid amplitude fluctuations). These effects are quantified by two critical observables:

- **Scintillation Bandwidth ($\Delta \nu$):** Characterizes rapid fluctuations in amplitude due to ISM turbulence.
- **Decorrelation Bandwidth ($\Delta f$):** Frequency range over which the signal's phase loses coherence.

In this project, we study FRB repeater signals from the M81 galaxy group, focusing on numerical computation of scintillation and decorrelation bandwidths using CHIME/Pulsar and CHIME/FRB systems. We also model potential quantum gravity corrections to spacetime geometry affecting FRB propagation.

---

## Theoretical Background

FRB signals propagating through turbulent ISM experience stochastic scattering effects. These effects are captured via the autocorrelation function (ACF), given by:

$$
\text{ACF}(\nu, \tau) = \int d\nu' \mathcal{I}(\nu + \nu') \mathcal{I}(\nu' + \tau),
$$

where:

- $\mathcal{I}$: Intensity of the radio signal.
- $\nu$: Frequency offset.
- $\tau$: Time lag.

The ACF's characteristic width along the frequency axis provides the scintillation bandwidth ($\Delta \nu$), while decorrelation bandwidth ($\Delta f$) is inferred from its rapid phase loss:

$$
\Delta \nu \sim \frac{1}{2\pi \tau_d}, \quad \Delta f \sim \frac{1}{2\pi D},
$$

where $\tau_d$ and $D$ describe scattering time and ISM turbulence, respectively.

FRBs may also offer insight into quantum gravity via spacetime corrections. For black hole-white hole tunneling, the tunneling amplitude modifies spacetime curvature, introducing frequency-dependent delays:

$$
\Delta t \sim \frac{\Delta E}{(1+z)^2 H_0},
$$

where $\Delta E$ is the energy range, $z$ is the redshift, and $H_0$ is the Hubble constant.

---

## Code Functionality

The primary code (*readScript_.py*) analyzes FRB-ISM interactions and quantum gravity implications using autocorrelation and spectral analysis:

1. **Autocorrelation Function (ACF) Computation:**
   - Numerically calculates ACFs from FRB intensity data.
   - Provides frequency-dependent scattering and scintillation bandwidths.

2. **Channel Flagging:**
   - Identifies and removes channels with amplitudes exceeding defined thresholds to reduce noise.

3. **Scintillation Bandwidth Estimation:**
   - Determines $\Delta \nu$ by fitting Gaussian profiles to ACF widths.

4. **Visualization:**
   - Plots overlapping signal densities, ACFs, and scintillation bandwidth fits.

5. **Quantum Gravity Corrections:**
   - Models FRB propagation using black-white hole tunneling amplitudes.

---

## Results and Visualizations

![Autocorrelation and Fit](https://github.com/IsolatedSingularity/FRB-Analysis/blob/main/Plots/ACF_Sum%3B%20Fit.png?raw=true)
*Autocorrelation function of an FRB signal with Gaussian fit for scintillation bandwidth estimation.*

The scintillation bandwidth varies across observed frequencies due to ISM turbulence. The ACFs reveal rapid amplitude fluctuations, enabling precise $\Delta \nu$ and $\Delta f$ calculations.

---

## Caveats

- **Extreme Limits of ACF:**
  - Near singularities, the ACF's behavior introduces degeneracies, leading to redundant computation. Optimizing these cases would improve computational efficiency.

- **Signal Quality and Noise:**
  - Channel flagging is sensitive to threshold values, potentially excluding valid data or retaining noise.

- **Simplified Quantum Gravity Models:**
  - Black-white hole tunneling models assume idealized geometry, omitting higher-dimensional corrections or exotic compact objects.

---

## Next Steps

- [x] Extend ACF analysis to include time-dependent scattering effects in the ISM.
- [ ] Explore multi-dimensional black-white wormhole tunneling models.
- [ ] Integrate machine learning for improved channel flagging and noise suppression.
- [ ] Test code scalability on larger FRB datasets using Canada Compute - Cedar.

---

> **Tip:** Always inspect flagged channels to ensure valid signal data is not excluded.

> **Note:** Detailed derivations of ACF formulas, ISM turbulence modeling, and quantum gravity corrections can be found in the associated PDF in the main repository.

