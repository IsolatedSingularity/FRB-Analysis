# Fast Radio Burst Cosmology
###### Under the supervision of [Professor Victoria Kaspi](https://www.physics.mcgill.ca/~vkaspi/) at [McGill University](https://www.mcgill.ca/). In collaboration with the [CHIME experiment](https://chime-experiment.ca/en) and [Canada Compute - Cedar](https://docs.alliancecan.ca/wiki/Cedar).

![alt text](https://github.com/IsolatedSingularity/FRB-Analysis/blob/main/Plots/ACF_Sum%3B%20Fit.png)

## Objective

In quantum gravity it is possible to describe a black hole quantum tunneling into a while hole in the classical limit [[1]](https://inspirehep.net/literature/1403709). This could occur in events where a black hole collapses faster than the Hawking radiation time, and so the matter in the black hole interior must tunnel to corresponding white hole interior. These tunnelling events are thought to be connected to with periodic radio signals occuring at low time scales, known as **fast radio bursts** (FRB) [[2]](https://www.researchgate.net/publication/265730144_Fast_Radio_Bursts_and_White_Hole_Signals). These signals are usually seen to be emitted by pulsars or merger systems and as the signal travels through the interstellar medium their signals experience dispersion causing different frequencies to arrive at different times. Moreover each are subject to spacetime dependent background noise given by primordial perturbations in the very early universe. This makes it interesting to study correlations of signals at different spatial separations that give us on how two sources of a signal can interact with one another at different distances. This is determined by statistical osbverables such as **decorrelation and scintillation bandwidths** which respectively quantify the range of frequencies over which the phases of a radio wave emitted by a source lose correlation during propagation and the rapid fluctuations in the amplitude of the radio signal during so.

We characterize this information for fast radio burst repeater signals detected from the direction of the M81 galaxy group. Numerically computed scintillation and decorrelation bandwidths of signals with the CHIME/Pulsar and CHIME/FRB systems amongst non-linear cosmological noise, after performing functional time series. We also analytically look at modelling this with with black-white hole tunneling events for corrections to the underlying spacetime geometry. 


## Code Functionality

*readScript_.py*  Fast Radio Burst (FRB) signals, particularly in the context of their interaction with the interstellar medium and their potential connection to quantum gravity phenomena. FRBs are brief bursts of radio waves from celestial sources, and this code helps researchers explore their properties.

The code begins by importing essential Python libraries and then reads filterbank data, which contains critical information about the frequency, time, and amplitude of the FRB signals. To enhance data quality, it identifies and removes channels with amplitudes exceeding defined thresholds.

The core of the code lies in its Auto Correlation Functions (ACF) calculations. It offers multiple methods to compute ACFs, each tailored to different analysis requirements. These ACFs provide insights into how FRB signals evolve as they propagate through spacetime. The code also fits the ACFs using mathematical models, notably Lorentzian functions, to extract valuable information.

Additionally, the code estimates the scintillation width, a key parameter characterizing the dynamic behavior of FRB signals in different frequency bands. It generates numerous plots to visualize the results, aiding in the comprehensive understanding of the data.

This analysis is particularly intriguing as it connects the study of FRBs to quantum gravity phenomena, where black holes may transition into white holes. By examining the interaction of FRB signals with their environment, researchers gain deeper insights into the mysteries of the cosmos, enhancing our understanding of astrophysical and quantum gravity phenomena.

The report covering all the theory and code can be found in the main repository as a PDF file.

## Caveats

The method in which points are detected within the wake is done using complex convex hulls. This algorithm
becomes problematic when the blown up deficit angle is replaced by its actual value of $\alpha = 8 \pi G \mu$ which
is very small and thus the wake becomes a plane. The algorithm is based on connecting simplices along
different vertices and does not work when the topology of the object is in 1D. Next, when converting from
physical to comoving coordinates, one uses an inverse scaling factor of the form $a^{−1}(z) = (1 − z)/z$, which
can also be substituted for $a^{−1}(t_0) \sim 10^3$ for current observations. This scaling becomes an issue when
wanting to scale physical axes to redshift axes using the numerical function from the astropy package, which
doesn’t converge for small $\mathcal{O}(1)$ or large $\mathcal{O}(1000)$ values of redshift. Thus, we are left with to work in a
snapshot of physical coordinates to substitute for a continuous comoving coordinate system.

## Next Steps

At this point, the problem statement has been defined and the foundational code base has been laid. Potential additions to the code would include higher dimensional topological defect signals, and including an algorithm to invert the redshift function without the limitation of convergence.
