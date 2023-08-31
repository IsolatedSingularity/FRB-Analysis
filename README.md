# FRB-Analysis
###### Under the supervision of [Professor Victoria Kaspi](https://www.physics.mcgill.ca/~vkaspi/) at [McGill University](https://www.mcgill.ca/). In collaboration with the [CHIME experiment](https://chime-experiment.ca/en) and [Canada Compute - Cedar](https://docs.alliancecan.ca/wiki/Cedar).

![alt text](https://github.com/IsolatedSingularity/FRB-Analysis/blob/main/Plots/ACF_Sum%3B%20Fit.png)

## Objective

In quantum gravity it is possible to describe a black hole (which formed from collapse) quantum tunneling into a while hole without violating classical general relativity [[1]](https://inspirehep.net/literature/1403709)

Characterized the statistics of fast radio burst repeater signals detected from the direction of the M81 galaxy group. Numerically computed scintillation and decorrelation bandwidths of signals with the CHIME/Pulsar and CHIME/FRB systems amongst non-linear cosmological noise.

The strings occur in a class of renormalizable quantum field theories and are stable under specific conditions on the spacetime manifold's homotopy groups. Their dynamics are given by the Nambu-Gotto action, much like bosonic excitations in string theory. What's more is that gravitational backreaction during the early universe causes primordial ΛCDM noise which hides the string signal. Thus the purpose of this repository is to develop statistics to efficiently extract the cosmic string signal admist the non-linear noise through the framework of 21cm cosmology.

## Code Functionality

*Cosmic String Extraction Statistics.py* builds the cosmic string signal from scratch as a finite density of energy radiating a certain temperature difference admist the 21cm background temperature map. This propagates through spacetime and traces out a wake which has the temperature gradient defined on its convex hull. Then using universe simulations done with 21cmFAST, the string signal is embedded in the primordial noise. Finally, to extract the dynamic signal we make use of statistics such as correlation functions, matched filters, and wavelets. The output are plots of these statistics when the signal of the string is detected in the noise.

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
