---
layout: page
permalink: /examples/
title: Demonstration examples for SReachTools
---

We consider various stochastic LTI systems like

* the relative dynamics of a spacecraft in space to a docking station
* the stochastic double integrator.

All these examples listed below were written using MATLAB's LiveEditor.
They are available in the `examples` folder of the project. 

**NOTE**: The HTML pages will open in the same page if CTRL/CMD key is not pressed.

1. **Forward stochastic reachability**: Computes the probability of the state of the spacecraft lying in a set at future time of interest with/without safety requirements. We also demonstrate how SReachTools can handle a deterministic and a stochastic initial state. [[PDF](https://github.com/abyvinod/SReachTools/raw/master/examples/forwardStochasticReachCWH.pdf)] [[HTML](https://htmlpreview.github.io/?https://github.com/abyvinod/SReachTools/blob/master/examples/forwardStochasticReachCWH.html)]
1. **Underapproximative verification using Fourier transforms and convex optimization**: Computes a lower bound on the safety probability for a spacecraft rendezvous and docking probability starting from a given state, its associated controller, and a subset of all the initial states from which a desired safety probability can be met. 
[[PDF](https://github.com/abyvinod/SReachTools/raw/master/examples/FtCVXUnderapproxVerifyCWH.pdf)] [[HTML](https://htmlpreview.github.io/?https://github.com/abyvinod/SReachTools/blob/master/examples/FtCVXUnderapproxVerifyCWH.html)]
1. **Stochastic reachability of a target tube using dynamic programming**: Computes a grid-based dynamic programming solution to the stochastic reachability of various target tubes for a stochastic double integrator. [[PDF](https://github.com/abyvinod/SReachTools/raw/master/examples/doubleIntegratorDynamicProgramming.pdf)] [[HTML](https://htmlpreview.github.io/?https://github.com/abyvinod/SReachTools/blob/master/examples/doubleIntegratorDynamicProgramming.html)]
1. **Underapproximative verification using Lagrangian (set operations) approach**: Coming soon

Please feel free to add requests for more examples at [https://github.com/abyvinod/SReachTools/issues](https://github.com/abyvinod/SReachTools/issues).

## Acknowledgements

Thanks to [https://htmlpreview.github.io/](https://htmlpreview.github.io/) for providing a way to link URL from GitHub.
