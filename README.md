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

```python
# Importing essential libraries
import numpy as np
from scipy import signal
import h5py

# Function to compute decorrelation bandwidth
def dedisperse(spectraData, dm=87.757, freq_hi=800.0, tsamp=0.00004096, foff=0.390625):
    spectraData_dedispersed = np.asarray(spectraData)
    freq_lo = freq_hi
    for channelIndex in range(np.shape(spectraData_dedispersed)[0]):
        tau = dm * (1/freq_lo**2 - 1/freq_hi**2)
        num_samples = int(round(-tau / tsamp))
        spectraData_dedispersed[channelIndex, :] = np.roll(spectraData_dedispersed[channelIndex, :], num_samples, axis=0)
        freq_lo -= foff
    return spectraData_dedispersed
```

1. **Autocorrelation Function (ACF) Computation:**
   - Numerically calculates ACFs from FRB intensity data.
   - Provides frequency-dependent scattering and scintillation bandwidths.

```python
# Define autocorrelation function (ACF)
def autocorrelation(data):
    result = signal.correlate(data, data, mode='full')
    return result[result.size // 2:]
```

2. **Channel Flagging:**
   - Identifies and removes channels with amplitudes exceeding defined thresholds to reduce noise.

```python
# Flagging noisy channels
def flag_channels(freq_spectrum, threshold=3):
    mean, std_dev = np.mean(freq_spectrum), np.std(freq_spectrum)
    flagged = [i for i, val in enumerate(freq_spectrum) if abs(val - mean) > threshold * std_dev]
    return flagged
```

3. **Scintillation Bandwidth Estimation:**
   - Determines $\Delta \nu$ by fitting Gaussian profiles to ACF widths.

```python
# Fit Gaussian to ACF width
def fit_gaussian(acf_data):
    popt, _ = signal.curve_fit(lambda x, a, b, c: a * np.exp(-(x-b)**2 / (2*c**2)),
                               np.arange(len(acf_data)), acf_data)
    return popt[1], popt[2]  # Mean and standard deviation
```

4. **Visualization:**
   - Plots overlapping signal densities, ACFs, and scintillation bandwidth fits.

```python
import matplotlib.pyplot as plt

# Visualization function
def plot_acf(acf_data):
    plt.plot(acf_data)
    plt.xlabel('Frequency Lag')
    plt.ylabel('Autocorrelation')
    plt.title('ACF of FRB Signal')
    plt.show()
```

5. **Quantum Gravity Corrections:**
   - Models FRB propagation using black-white hole tunneling amplitudes.

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

> [!TIP]
> Always inspect flagged channels to ensure valid signal data is not excluded.

> [!NOTE]
> Detailed derivations of ACF formulas, ISM turbulence modeling, and quantum gravity corrections can be found in the associated PDF in the main repository.
