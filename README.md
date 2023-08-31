# Fast Radio Burst Tunneling
###### Under the supervision of [Professor Victoria Kaspi](https://www.physics.mcgill.ca/~vkaspi/) at [McGill University](https://www.mcgill.ca/). In collaboration with the [CHIME experiment](https://chime-experiment.ca/en) and [Canada Compute - Cedar](https://docs.alliancecan.ca/wiki/Cedar).

![alt text](https://github.com/IsolatedSingularity/FRB-Analysis/blob/main/Plots/ACF_Sum%3B%20Fit.png)

## Objective

In quantum gravity it is possible to describe a black hole quantum tunneling into a while hole in the classical limit [[1]](https://inspirehep.net/literature/1403709). This could occur in events where a black hole collapses faster than the Hawking radiation time, and so the matter in the black hole interior must tunnel to corresponding white hole interior. These tunnelling events are thought to be connected to with periodic radio signals occuring at low time scales, known as **fast radio bursts** (FRB) [[2]](https://www.researchgate.net/publication/265730144_Fast_Radio_Bursts_and_White_Hole_Signals). These signals are usually seen to be emitted by pulsars or merger systems and as the signal travels through the interstellar medium their signals experience dispersion causing different frequencies to arrive at different times. Moreover each are subject to spacetime dependent background noise given by primordial perturbations in the very early universe. This makes it interesting to study correlations of signals at different spatial separations as it give us information on how two sources of a signal can interact with one another. This is determined by statistical osbverables such as **decorrelation and scintillation bandwidths** which respectively quantify the range of frequencies over which the phases of a radio wave emitted by a source lose correlation during propagation and the rapid fluctuations in the amplitude of the radio signal during so.

We characterize this information for fast radio burst (FRB) repeater signals detected from the direction of the M81 galaxy group. Numerically computed scintillation and decorrelation bandwidths of signals with the CHIME/Pulsar and CHIME/FRB systems amongst non-linear cosmological noise, after performing functional time series. We also analytically look at modelling this with with black-white hole tunneling events for corrections to the underlying spacetime geometry. 


## Code Functionality

*readScript_.py* characterizes interaction between the fast radio bursts their their interstellar medium and their potential connection to quantum gravity phenomena. The core of the code lies in its auto-correlated functionals (ACFs) of frequencies at different points in spacetime. It offers multiple methods to compute ACFs, each tailored to different analysis requirements. These ACFs provide insights into how FRB signals evolve as they propagate through spacetime. To enhance data quality, it identifies and removes channels with amplitudes exceeding defined thresholds. The code thereafter estimates the decorrelation bandwidth, a key parameter characterizing the dynamic behavior of FRB signals in different frequency bands. It generates numerous plots to visualize the results in the form of interacting signal densities and how the scintillation bandwidth can be computed from the overlapping densities.

## Caveats

One caveat to consider is the potential behavior of the frequency functional at extreme limits. In cases where the functional approaches singularities, the code has degeneracies in the information bordering the signal distribution. Removing these degenracies to avoid redudnant computation is essential for scalability in terms of its time complexity.

## Next Steps

At this point, the problem statement has been defined and the foundational code base has been laid. Potential additions to the code would include higher dimensional black-white wormhole tunneling events with sources that are causally disconnected. It would also be informative to study this in the context of a foliated spacetime, such as in the ADM formalism.
