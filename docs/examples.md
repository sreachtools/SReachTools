---
layout: page
permalink: /examples/
title: Demonstration examples for SReachTools
---

In these set of examples, we apply SReachTools to perform forward (prediction of the stochasticity of the state in future) and backward (verification and controller synthesis) stochastic reachability in various stochastic LTI systems like:

* the relative dynamics of a spacecraft in space to a docking station,
* the stochastic double integrator, and
* the automated anesthesia delivery system.

Please feel free to add requests for more examples at [https://github.com/unm-hscl/SReachTools/issues](https://github.com/unm-hscl/SReachTools/issues).
We used MATLAB's LiveEditor to write the following examples.
They are available in the `examples` folder of the project. 

{% include important-note.html content="HTML examples are not mobile friendly. We're working on it." %}

1. **Forward stochastic reachability**: Computes the probability of the state of the spacecraft lying in a set at future time of interest with/without safety requirements. We also demonstrate how SReachTools can handle a deterministic and a stochastic initial state. [[PDF](https://github.com/unm-hscl/SReachTools/raw/master/examples/forwardStochasticReachCWH.pdf)] [[HTML](forwardStochasticReachCWH.html)] 
1. **Underapproximative verification using Fourier transforms and convex optimization**: 
    1. **Spacecraft rendezvous and docking problem**: Computes a lower bound on the safety probability for a spacecraft rendezvous and docking probability starting from a given state, its associated controller, and a subset of all the initial states from which a desired safety probability can be met.
[[PDF](https://github.com/unm-hscl/SReachTools/raw/master/examples/FtCVXUnderapproxVerifyCWH.pdf)] [[HTML](FtCVXUnderapproxVerifyCWH.html)]
    1. **Automated anesthesia delivery system**: Computes a subset of all the initial states from which a depth of hypnosis remains in a desired range (probalistically) by the automated anesthesia delivery system under process uncertainty.
[[PDF](https://github.com/unm-hscl/SReachTools/raw/master/examples/AutomatedAnesthesiaDelivery.pdf)] [[HTML](AutomatedAnesthesiaDelivery.html)]
1. **Stochastic reachability of a target tube using dynamic programming**: Computes a grid-based dynamic programming solution to the stochastic reachability of various target tubes for a stochastic double integrator. [[PDF](https://github.com/unm-hscl/SReachTools/raw/master/examples/doubleIntegratorDynamicProgramming.pdf)] [[HTML](doubleIntegratorDynamicProgramming.html)]
1. **Underapproximative verification using Lagrangian (set operations) approach**: Computes an underapproximation of the stochastic reach-avoid set using Lagrangian, set-based, methods. Grid free recursion process implemented with the
Model Parametric Toolbox (MPT). [[PDF](https://github.com/unm-hscl/SReachTools/raw/master/examples/lagrangianApproximations.pdf)] [[HTML](lagrangianApproximations.html)]

<!-- Add {:target="_blank"} if it is desired that the page opens in a new window.-->
