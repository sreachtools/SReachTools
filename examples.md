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

{% include important-note.html content="HTML examples are not mobile friendly. We're working on it." %}

1. **Forward stochastic reachability**: Computes the probability of the state of the spacecraft lying in a set at future time of interest with/without safety requirements. We also demonstrate how SReachTools can handle a deterministic and a stochastic initial state. [[PDF](https://github.com/abyvinod/SReachTools/raw/master/examples/forwardStochasticReachCWH.pdf)] [[HTML](forwardStochasticReachCWH.html)] 
1. **Underapproximative verification using Fourier transforms and convex optimization**: Computes a lower bound on the safety probability for a spacecraft rendezvous and docking probability starting from a given state, its associated controller, and a subset of all the initial states from which a desired safety probability can be met.
[[PDF](https://github.com/abyvinod/SReachTools/raw/master/examples/FtCVXUnderapproxVerifyCWH.pdf)] [[HTML](FtCVXUnderapproxVerifyCWH.html)]
1. **Stochastic reachability of a target tube using dynamic programming**: Computes a grid-based dynamic programming solution to the stochastic reachability of various target tubes for a stochastic double integrator. [[PDF](https://github.com/abyvinod/SReachTools/raw/master/examples/doubleIntegratorDynamicProgramming.pdf)] [[HTML](doubleIntegratorDynamicProgramming.html)]
1. **Underapproximative verification using Lagrangian (set operations) approach**: Computes an underapproximation of the stochastic reach-avoid set using Lagrangian, set-based, methods. Grid free recursion process implemented with the
Model Parametric Toolbox (MPT). [[PDF](https://github.com/abyvinod/SReachTools/raw/master/examples/lagrangianApproximations.pdf)] [[HTML](lagrangianApproximations.html)]

Please feel free to add requests for more examples at [https://github.com/abyvinod/SReachTools/issues](https://github.com/abyvinod/SReachTools/issues).

<!-- Add {:target="_blank"} if it is desired that the page opens in a new window.-->
